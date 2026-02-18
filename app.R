# app.R
# Incucyte Multi-File Import + Channel Normalization (Passage as surrogate BioRep)
# + Condition header parser -> receptor + treatment factors
#
# Features:
# - Imports Incucyte "wide-by-condition" exports (header line: "Date Time<TAB>Elapsed...")
# - Per-file channel assignment; default: if filename contains "red" anywhere -> NIR
# - Passage parsed from file header and stored as factor (surrogate biological replicate)
# - Normalization within Passage (ratio or log2 ratio; optional baseline normalize)
# - Adds receptor + treatment factors guessed from condition headers
# - Plot tab left of Join checks
# - Downloads: tidy long, joined wide, normalized, stats-ready long

library(shiny)
library(tidyverse)
library(readr)
library(stringr)

# ---------------------------
# Helpers: parse Incucyte header metadata
# ---------------------------
find_data_header_row <- function(lines) {
  idx <- which(str_detect(lines, "^Date Time\\tElapsed\\b"))
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

extract_key_value <- function(lines, key) {
  pat <- paste0("^", stringr::fixed(key), "\\s*:\\s*(.*)$")
  hit <- lines[str_detect(lines, pat)]
  if (length(hit) == 0) return(NA_character_)
  val <- str_match(hit[1], pat)[, 2]
  if (is.na(val)) NA_character_ else str_trim(val)
}

read_incucyte_header_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  header_i <- find_data_header_row(lines)
  if (is.na(header_i)) stop("Couldn't find data header row starting with 'Date Time<TAB>Elapsed'.")
  
  meta_lines <- lines[1:(header_i - 1)]
  tibble(
    vessel_name = extract_key_value(meta_lines, "Vessel Name"),
    metric      = extract_key_value(meta_lines, "Metric"),
    cell_type   = extract_key_value(meta_lines, "Cell Type"),
    passage     = extract_key_value(meta_lines, "Passage"),
    notes       = extract_key_value(meta_lines, "Notes"),
    analysis    = extract_key_value(meta_lines, "Analysis")
  )
}

# ---------------------------
# Reader: Incucyte wide-by-condition export -> tidy long
# ---------------------------
read_incucyte_wide_conditions <- function(path, drop_stderr = TRUE) {
  lines <- readLines(path, warn = FALSE)
  
  header_i <- find_data_header_row(lines)
  if (is.na(header_i)) stop("Couldn't find header row starting with 'Date Time<TAB>Elapsed'.")
  
  txt <- paste(lines[header_i:length(lines)], collapse = "\n")
  
  dat <- read_delim(
    I(txt),
    delim = "\t",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )
  
  if (ncol(dat) < 3) stop("Parsed <3 columns. Is this the right Incucyte export format?")
  names(dat)[1:2] <- c("datetime", "elapsed")
  
  if (drop_stderr) dat <- dat %>% select(-matches("\\(Std Err"))
  
  dat %>%
    mutate(
      datetime = as.character(datetime),
      elapsed  = suppressWarnings(as.numeric(elapsed))
    ) %>%
    pivot_longer(
      cols = -c(datetime, elapsed),
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(
      condition = as.character(condition),
      value     = suppressWarnings(as.numeric(value))
    )
}

# ---------------------------
# Condition header parsing: receptor + treatment factors
# ---------------------------

canonicalize_receptor_token <- function(x) {
  z <- x %>%
    str_to_lower() %>%
    str_replace_all("[α]", "a") %>%
    str_replace_all("[β]", "b") %>%
    str_replace_all("[^a-z0-9\\s]", " ") %>%
    str_squish()
  
  dplyr::case_when(
    z %in% c("era", "er a", "esr1", "er1") ~ "era",
    z %in% c("erb", "er b", "esr2", "er2") ~ "erb",
    z %in% c("pgr", "pr", "prg") ~ "pgr",
    z %in% c("pra", "pr a", "pr a isoform", "pr a iso", "pr-a", "pr a iso") ~ "pra",
    z %in% c("prb", "pr b", "pr b isoform", "pr b iso", "pr-b", "pr b iso") ~ "prb",
    z %in% c("gr", "nr3c1", "glucocorticoid receptor") ~ "gr",
    TRUE ~ NA_character_
  )
}

is_vehicle_token <- function(x) {
  z <- str_to_lower(str_squish(x))
  z %in% c("veh", "vehicle", "etoh", "ethanol", "dmso", "pbs", "media", "control")
}

parse_treatment_token <- function(x) {
  z <- str_squish(x)
  
  m <- str_match(
    z,
    regex("^\\s*([0-9]+\\.?[0-9]*)\\s*(pm|nm|um|µm|mm)\\s*([a-z0-9\\-]+)\\s*$",
          ignore_case = TRUE)
  )
  if (all(is.na(m))) return(NA_character_)
  
  amount <- m[, 2]
  unit   <- toupper(m[, 3])
  unit   <- ifelse(unit == "µM", "uM", unit)
  comp   <- toupper(m[, 4])
  
  paste0(amount, unit, " ", comp)
}

canonicalize_receptor_combo <- function(receptors) {
  order <- c("era", "erb", "pgr", "pra", "prb", "gr")
  recs <- unique(na.omit(receptors))
  if (length(recs) == 0) return(NA_character_)
  recs <- recs[order(match(recs, order, nomatch = 999))]
  paste(recs, collapse = " + ")
}

guess_receptor_treatment <- function(condition_vec) {
  tibble(condition_raw = condition_vec) %>%
    mutate(
      condition_norm = condition_raw %>%
        str_replace_all("\\s*\\+\\s*", " + ") %>%
        str_squish(),
      tokens = str_split(condition_norm, " \\+ ")
    ) %>%
    mutate(
      parsed = purrr::map(tokens, function(tok) {
        tok <- str_squish(tok)
        
        receptor_tokens <- purrr::map_chr(tok, canonicalize_receptor_token)
        receptors <- receptor_tokens[!is.na(receptor_tokens)]
        non_receptor <- tok[is.na(receptor_tokens)]
        
        veh_idx <- purrr::map_lgl(non_receptor, is_vehicle_token)
        conc_parsed <- purrr::map_chr(non_receptor, parse_treatment_token)
        conc_idx <- !is.na(conc_parsed)
        
        treatment <- dplyr::case_when(
          any(conc_idx) ~ paste(unique(conc_parsed[conc_idx]), collapse = " + "),
          any(veh_idx) ~ "VEH",
          length(non_receptor) == 0 ~ NA_character_,
          TRUE ~ paste(non_receptor, collapse = " + ")
        )
        
        list(
          receptor = canonicalize_receptor_combo(receptors),
          treatment = treatment
        )
      })
    ) %>%
    tidyr::unnest_wider(parsed) %>%
    mutate(
      receptor = factor(receptor),
      treatment = factor(treatment)
    ) %>%
    select(condition_raw, receptor, treatment)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Incucyte Multi-File Import + Channel Normalization (Passage + receptor/treatment parsing)"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "files",
        "Upload Incucyte export files (.txt/.tsv/.csv)",
        multiple = TRUE,
        accept = c(".txt", ".tsv", ".csv")
      ),
      checkboxInput("drop_stderr", "Drop '(Std Err ...)' columns", value = TRUE),
      uiOutput("channel_map_ui"),
      hr(),
      uiOutput("norm_ui"),
      actionButton("run", "Import + Process", class = "btn-primary"),
      hr(),
      downloadButton("download_long", "Download: tidy long (csv)"),
      downloadButton("download_wide", "Download: joined wide (csv)"),
      downloadButton("download_norm", "Download: normalized (csv)"),
      downloadButton("download_stats", "Download: stats-ready long (csv)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot", height = 420)),
        tabPanel("Join checks", verbatimTextOutput("join_checks")),
        tabPanel("Raw (tidy long)", tableOutput("preview_long")),
        tabPanel("Joined (wide by channel)", tableOutput("preview_wide")),
        tabPanel("Normalized", tableOutput("preview_norm")),
        tabPanel("Stats-ready", tableOutput("preview_stats"))
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  
  # Dynamic UI: per-file channel assignment
  output$channel_map_ui <- renderUI({
    req(input$files)
    fns <- input$files$name
    
    tagList(
      h4("Assign a fluorescence channel to each file"),
      tags$p(tags$small("Default: if filename contains 'red' anywhere → NIR.")),
      lapply(seq_along(fns), function(i) {
        fname <- fns[i]
        
        default <- if (str_detect(tolower(fname), "red")) "NIR" else
          if (str_detect(tolower(fname), "nir")) "NIR" else
            if (str_detect(tolower(fname), "gfp|green")) "GFP" else
              if (str_detect(tolower(fname), "orange")) "Orange" else
                "Other"
        
        fluidRow(
          column(8, tags$small(fname)),
          column(
            4,
            selectInput(
              inputId = paste0("chan_", i),
              label = NULL,
              choices = c("GFP", "NIR", "Orange", "Red", "Other"),
              selected = default
            )
          )
        )
      })
    )
  })
  
  # Safe getter for dynamic inputs
  channel_map <- reactive({
    req(input$files)
    fns <- input$files$name
    
    get_chan <- function(i, fname) {
      id <- paste0("chan_", i)
      val <- input[[id]]
      if (is.null(val) || length(val) == 0 || is.na(val) || val == "") {
        if (str_detect(tolower(fname), "red")) return("NIR")
        if (str_detect(tolower(fname), "nir")) return("NIR")
        if (str_detect(tolower(fname), "gfp|green")) return("GFP")
        if (str_detect(tolower(fname), "orange")) return("Orange")
        return("Other")
      }
      val
    }
    
    tibble(
      file = fns,
      path = input$files$datapath,
      channel = purrr::map_chr(seq_along(fns), ~ get_chan(.x, fns[.x]))
    )
  })
  
  # Normalization UI
  output$norm_ui <- renderUI({
    req(channel_map())
    chans <- sort(unique(channel_map()$channel))
    if (length(chans) == 0) return(NULL)
    
    tagList(
      h4("Normalization"),
      selectInput("signal_channel", "Signal channel", choices = chans, selected = chans[1]),
      selectInput(
        "control_channel",
        "Control channel (e.g., NIR for transfection efficiency)",
        choices = chans,
        selected = if ("NIR" %in% chans) "NIR" else chans[min(2, length(chans))]
      ),
      radioButtons(
        "norm_method",
        "Method",
        choices = c(
          "None (keep raw signal channel)" = "none",
          "Ratio (signal/control)" = "ratio",
          "Log2 ratio (log2(signal/control))" = "log2ratio"
        ),
        selected = "ratio"
      ),
      checkboxInput(
        "baseline_norm",
        "Also baseline-normalize within Passage+Condition (divide by first non-NA timepoint)",
        value = FALSE
      )
    )
  })
  
  # Import all files to tidy long + parse header metadata (incl Passage) + add receptor/treatment
  raw_long <- eventReactive(input$run, {
    cm <- channel_map()
    
    dat <- purrr::pmap_dfr(cm, function(file, path, channel) {
      meta <- read_incucyte_header_meta(path)
      dat0 <- read_incucyte_wide_conditions(path, drop_stderr = isTRUE(input$drop_stderr))
      
      passage_chr <- meta$passage[[1]]
      passage_lab <- if (!is.na(passage_chr) && passage_chr != "") paste0("Passage_", passage_chr) else "Passage_NA"
      
      dat0 %>%
        mutate(
          file = file,
          channel = channel,
          passage = factor(passage_lab),
          vessel_name = meta$vessel_name[[1]] %||% NA_character_,
          metric      = meta$metric[[1]] %||% NA_character_,
          cell_type   = meta$cell_type[[1]] %||% NA_character_,
          analysis    = meta$analysis[[1]] %||% NA_character_
        )
    }) %>%
      select(file, channel, passage, vessel_name, metric, cell_type, analysis,
             datetime, elapsed, condition, value)
    
    # Add receptor/treatment factors based on condition headers
    annot <- guess_receptor_treatment(unique(dat$condition))
    dat %>%
      left_join(annot, by = c("condition" = "condition_raw")) %>%
      relocate(receptor, treatment, .after = condition)
  }, ignoreInit = TRUE)
  
  output$preview_long <- renderTable({
    req(raw_long())
    head(raw_long(), 20)
  })
  
  # Join checks
  output$join_checks <- renderPrint({
    req(raw_long())
    rl <- raw_long()
    
    files_loaded <- rl %>%
      distinct(file, channel, passage, vessel_name, metric, cell_type, analysis) %>%
      arrange(passage, channel, file)
    
    per_channel <- rl %>%
      group_by(channel) %>%
      summarize(
        n_rows = n(),
        n_passages = n_distinct(passage),
        n_conditions = n_distinct(condition),
        n_times = n_distinct(elapsed),
        .groups = "drop"
      )
    
    # parsing coverage
    parse_cov <- rl %>%
      summarize(
        n_conditions = n_distinct(condition),
        n_receptor_missing = sum(is.na(receptor)),
        n_treatment_missing = sum(is.na(treatment))
      )
    
    list(
      files_loaded = files_loaded,
      per_channel_summary = per_channel,
      parsing_coverage = parse_cov,
      note = "Passage parsed from file header line 'Passage: <n>' and stored as factor surrogate biological replicate."
    )
  })
  
  # Wide join WITHIN Passage (stats-friendly)
  wide_joined_passage <- reactive({
    req(raw_long())
    raw_long() %>%
      select(passage, channel, elapsed, condition, receptor, treatment, value) %>%
      group_by(passage, channel, elapsed, condition, receptor, treatment) %>%
      summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = channel, values_from = value)
  })
  
  output$preview_wide <- renderTable({
    req(wide_joined_passage())
    head(wide_joined_passage(), 20)
  })
  
  # Normalized within Passage
  normalized_passage <- reactive({
    req(wide_joined_passage())
    w <- wide_joined_passage()
    
    method <- input$norm_method %||% "ratio"
    sig <- input$signal_channel
    ctl <- input$control_channel
    
    if (is.null(sig) || is.null(ctl) || !(sig %in% names(w)) || !(ctl %in% names(w))) {
      return(w %>% mutate(value_norm = NA_real_))
    }
    
    out <- w %>%
      mutate(
        value_norm = case_when(
          method == "none" ~ .data[[sig]],
          method == "ratio" ~ .data[[sig]] / .data[[ctl]],
          method == "log2ratio" ~ log2(.data[[sig]] / .data[[ctl]]),
          TRUE ~ NA_real_
        )
      )
    
    if (isTRUE(input$baseline_norm)) {
      out <- out %>%
        group_by(passage, condition) %>%
        mutate(
          baseline = value_norm[which(!is.na(value_norm))[1]],
          value_norm = value_norm / baseline
        ) %>%
        ungroup() %>%
        select(-baseline)
    }
    
    out
  })
  
  output$preview_norm <- renderTable({
    req(normalized_passage())
    head(normalized_passage(), 20)
  })
  
  # Stats-ready long table: (Passage surrogate replicate, receptor, treatment, condition, time)
  stats_long <- reactive({
    req(normalized_passage())
    normalized_passage() %>%
      select(passage, elapsed, condition, receptor, treatment, value_norm) %>%
      arrange(passage, receptor, treatment, condition, elapsed)
  })
  
  output$preview_stats <- renderTable({
    req(stats_long())
    head(stats_long(), 20)
  })
  
  # Plot: thin lines = each Passage trajectory; thick line = mean across Passage
  output$plot <- renderPlot({
    req(stats_long())
    df <- stats_long()
    
    if (all(is.na(df$value_norm))) {
      plot.new()
      text(0.5, 0.5, "No normalized values to plot yet.\nUpload files, assign channels, choose signal/control, click 'Import + Process'.")
      return()
    }
    
    mean_df <- df %>%
      group_by(elapsed, condition, receptor, treatment) %>%
      summarize(mean = mean(value_norm, na.rm = TRUE), .groups = "drop")
    
    ggplot(df, aes(x = elapsed, y = value_norm, group = interaction(passage, condition))) +
      geom_line(alpha = 0.35) +
      geom_line(data = mean_df, aes(y = mean, group = condition), linewidth = 1.0) +
      labs(
        x = "Elapsed time",
        y = "Value (normalized)",
        title = "Normalized trajectories by Passage (thin) + mean across Passage (thick)"
      ) +
      theme_minimal()
  })
  
  # ---------------------------
  # Downloads
  # ---------------------------
  output$download_long <- downloadHandler(
    filename = function() paste0("incucyte_tidy_long_", Sys.Date(), ".csv"),
    content = function(file) {
      req(raw_long())
      write_csv(raw_long(), file)
    }
  )
  
  output$download_wide <- downloadHandler(
    filename = function() paste0("incucyte_joined_wide_by_passage_", Sys.Date(), ".csv"),
    content = function(file) {
      req(wide_joined_passage())
      write_csv(wide_joined_passage(), file)
    }
  )
  
  output$download_norm <- downloadHandler(
    filename = function() paste0("incucyte_normalized_by_passage_", Sys.Date(), ".csv"),
    content = function(file) {
      req(normalized_passage())
      write_csv(normalized_passage(), file)
    }
  )
  
  output$download_stats <- downloadHandler(
    filename = function() paste0("incucyte_stats_ready_long_", Sys.Date(), ".csv"),
    content = function(file) {
      req(stats_long())
      write_csv(stats_long(), file)
    }
  )
}

shinyApp(ui, server)