# app.R
# Incucyte Multi-File Import + Channel Normalization
# - Works with "wide-by-condition" Incucyte exports like your exp97.1_greenintensity_data.txt
#   (metadata header lines; real header begins with: "Date Time<TAB>Elapsed<...conditions...>")
# - Per-file channel assignment (defaults: if filename contains "red" -> NIR)
# - Normalization (signal/control) computed *within Passage* (Passage parsed from file header)
# - Plot tab is left of Join checks
# - Stores Passage as factor as surrogate biological replicate for stats

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
  # matches e.g. "Passage: 13"
  pat <- paste0("^", stringr::fixed(key), "\\s*:\\s*(.*)$")
  hit <- lines[str_detect(lines, pat)]
  if (length(hit) == 0) return(NA_character_)
  val <- str_match(hit[1], pat)[,2]
  if (is.na(val)) NA_character_ else str_trim(val)
}

read_incucyte_header_meta <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  header_i <- find_data_header_row(lines)
  if (is.na(header_i)) {
    stop("Couldn't find data header row starting with 'Date Time<TAB>Elapsed'.")
  }
  
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
  if (is.na(header_i)) {
    stop("Couldn't find header row starting with 'Date Time<TAB>Elapsed'.")
  }
  
  txt <- paste(lines[header_i:length(lines)], collapse = "\n")
  
  dat <- read_delim(
    I(txt),
    delim = "\t",
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  )
  
  if (ncol(dat) < 3) stop("Parsed <3 columns. Is this the right Incucyte export format?")
  names(dat)[1:2] <- c("datetime", "elapsed")
  
  if (drop_stderr) {
    dat <- dat %>% select(-matches("\\(Std Err"))
  }
  
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
# UI
# ---------------------------

ui <- fluidPage(
  titlePanel("Incucyte Multi-File Import + Channel Normalization (Passage = BioRep surrogate)"),
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
        
        # IMPORTANT: check "red" first -> NIR (per your request)
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
        # fallback defaults (same logic)
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
      channel = map_chr(seq_along(fns), ~ get_chan(.x, fns[.x]))
    )
  })
  
  # Normalization UI (appears once channels are known)
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
  
  # Import all files to tidy long + header metadata (incl Passage)
  raw_long <- eventReactive(input$run, {
    cm <- channel_map()
    
    purrr::pmap_dfr(cm, function(file, path, channel) {
      meta <- read_incucyte_header_meta(path)  # passage, vessel_name, etc.
      dat  <- read_incucyte_wide_conditions(path, drop_stderr = isTRUE(input$drop_stderr))
      
      passage_chr <- meta$passage[[1]]
      # store as factor label like "Passage_13" (handles NA gracefully)
      passage_lab <- if (!is.na(passage_chr) && passage_chr != "") paste0("Passage_", passage_chr) else "Passage_NA"
      
      dat %>%
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
  }, ignoreInit = TRUE)
  
  output$preview_long <- renderTable({
    req(raw_long())
    head(raw_long(), 20)
  })
  
  # Join checks summary
  output$join_checks <- renderPrint({
    req(raw_long())
    rl <- raw_long()
    
    per_file <- rl %>%
      distinct(file, channel, passage, vessel_name, metric, cell_type, analysis)
    
    by_chan <- rl %>%
      group_by(channel) %>%
      summarize(
        n_rows = n(),
        n_passages = n_distinct(passage),
        n_conditions = n_distinct(condition),
        n_times = n_distinct(elapsed),
        min_time = suppressWarnings(min(elapsed, na.rm = TRUE)),
        max_time = suppressWarnings(max(elapsed, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # pairwise condition overlap between channels (global)
    chans <- sort(unique(rl$channel))
    conds_by <- split(unique(rl$condition), rl$channel)
    
    pairwise <- tibble()
    if (length(chans) >= 2) {
      pairwise <- purrr::map_dfr(combn(chans, 2, simplify = FALSE), function(pair) {
        a <- pair[1]; b <- pair[2]
        ia <- conds_by[[a]] %||% character()
        ib <- conds_by[[b]] %||% character()
        tibble(
          channel_a = a,
          channel_b = b,
          conditions_in_a_not_b = length(setdiff(ia, ib)),
          conditions_in_b_not_a = length(setdiff(ib, ia)),
          conditions_in_both = length(intersect(ia, ib))
        )
      })
    }
    
    list(
      files_loaded = per_file,
      per_channel_summary = by_chan,
      pairwise_condition_overlap = pairwise,
      note = "Passage is parsed from file header line 'Passage: <n>' and stored as a factor surrogate for biological replicate."
    )
  })
  
  # Wide join WITHIN Passage (key for stats)
  wide_joined_passage <- reactive({
    req(raw_long())
    raw_long() %>%
      select(passage, channel, elapsed, condition, value) %>%
      group_by(passage, channel, elapsed, condition) %>%
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
  
  # Stats-ready long table: one row per (Passage, Condition, Time)
  stats_long <- reactive({
    req(normalized_passage())
    normalized_passage() %>%
      select(passage, elapsed, condition, value_norm) %>%
      arrange(passage, condition, elapsed)
  })
  
  output$preview_stats <- renderTable({
    req(stats_long())
    head(stats_long(), 20)
  })
  
  # Plot:
  # - thin lines = each Passage trajectory (surrogate biological replicate)
  # - thick line = mean across Passage
  output$plot <- renderPlot({
    req(stats_long())
    df <- stats_long()
    
    if (all(is.na(df$value_norm))) {
      plot.new()
      text(0.5, 0.5, "No normalized values to plot yet.\nUpload files, assign channels, choose signal/control, click 'Import + Process'.")
      return()
    }
    
    mean_df <- df %>%
      group_by(elapsed, condition) %>%
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