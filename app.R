# app.R
# Incucyte Multi-File Import + Channel Normalization
# - Passage parsed from header as surrogate biological replicate
# - Auto-parse receptor + treatment from condition headers
# - Plate size auto-detected from number of unique condition headers
# - Factor editor shown in plate-equivalent layout using rhandsontable:
#     * Condition header plate (read-only)
#     * Receptor plate (editable)
#     * Treatment plate (editable)
# - Prism export with one row per elapsed time
# - Plot defaults: facet by Passage; color by receptor

library(shiny)
library(tidyverse)
library(readr)
library(stringr)
library(rhandsontable)

# ---------------------------
# Helpers
# ---------------------------
`%||%` <- function(x, y) if (is.null(x)) y else x

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
      elapsed = suppressWarnings(as.numeric(elapsed))
    ) %>%
    pivot_longer(
      cols = -c(datetime, elapsed),
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(
      condition = as.character(condition),
      value = suppressWarnings(as.numeric(value))
    )
}

# Strict parser for condition headers only
canonicalize_receptor_token <- function(x) {
  z <- x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[α]", "a") %>%
    str_replace_all("[β]", "b") %>%
    str_replace_all("[^a-z0-9\\s]", " ") %>%
    str_squish()
  
  dplyr::case_when(
    z %in% c("era", "er a", "esr1", "er1") ~ "era",
    z %in% c("erb", "er b", "esr2", "er2") ~ "erb",
    z %in% c("pgr", "pr", "prg") ~ "pgr",
    z %in% c("pra", "pr a", "pr a isoform", "pr a iso") ~ "pra",
    z %in% c("prb", "pr b", "pr b isoform", "pr b iso") ~ "prb",
    z %in% c("ar", "androgen receptor", "nr3c4") ~ "AR",
    z %in% c("gr", "glucocorticoid receptor", "nr3c1") ~ "GR",
    z %in% c("mr", "mineralocorticoid receptor", "nr3c2") ~ "MR",
    TRUE ~ NA_character_
  )
}

# Looser parser for manual edits in the factor table
canonicalize_receptor_edit <- function(x) {
  z <- x %>%
    str_trim() %>%
    str_to_lower() %>%
    str_replace_all("[α]", "a") %>%
    str_replace_all("[β]", "b") %>%
    str_replace_all("[^a-z0-9\\s]", " ") %>%
    str_squish()
  
  dplyr::case_when(
    z %in% c("", "none", "na") ~ "none",
    z %in% c("era", "er a", "esr1", "er1") ~ "era",
    z %in% c("erb", "er b", "esr2", "er2") ~ "erb",
    z %in% c("pgr", "pr", "prg") ~ "pgr",
    z %in% c("pra", "pr a", "pr a isoform", "pr a iso") ~ "pra",
    z %in% c("prb", "pr b", "pr b isoform", "pr b iso") ~ "prb",
    z %in% c("ar", "androgen receptor", "nr3c4") ~ "AR",
    z %in% c("gr", "glucocorticoid receptor", "nr3c1") ~ "GR",
    z %in% c("mr", "mineralocorticoid receptor", "nr3c2") ~ "MR",
    TRUE ~ str_trim(x)
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
  order <- c("none", "era", "erb", "pgr", "pra", "prb", "AR", "GR", "MR")
  recs <- unique(na.omit(receptors))
  if (length(recs) == 0) return("none")
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
          length(non_receptor) == 0 ~ "VEH",
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

preview_file_lines <- function(path, n = 10) {
  lines <- readLines(path, warn = FALSE, n = n)
  tibble(line_no = seq_along(lines), line = lines)
}

ols_control_adjust <- function(df, sig_col, ctl_col) {
  x <- df[[ctl_col]]
  y <- df[[sig_col]]
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3) {
    df$value_norm <- NA_real_
    return(df)
  }
  
  fit <- stats::lm(y[ok] ~ x[ok])
  b <- unname(stats::coef(fit)[2])
  xbar <- mean(x[ok], na.rm = TRUE)
  
  adj <- rep(NA_real_, length(y))
  adj[ok] <- y[ok] - b * (x[ok] - xbar)
  
  df$value_norm <- adj
  df
}

# ---------------------------
# Plate layout helpers
# ---------------------------
detect_plate_layout <- function(n) {
  standard <- tibble(
    plate_size = c(6, 12, 24, 48, 96),
    n_rows = c(2, 3, 4, 6, 8),
    n_cols = c(3, 4, 6, 8, 12)
  )
  
  hit <- standard %>% filter(plate_size >= n) %>% slice(1)
  
  if (nrow(hit) == 0) {
    n_rows <- ceiling(sqrt(n))
    n_cols <- ceiling(n / n_rows)
    tibble(
      plate_size = n,
      n_rows = n_rows,
      n_cols = n_cols
    )
  } else {
    hit
  }
}

make_plate_df <- function(values, n_rows, n_cols, fill = "") {
  total <- n_rows * n_cols
  vals <- c(as.character(values), rep(fill, max(0, total - length(values))))
  mat <- matrix(vals, nrow = n_rows, ncol = n_cols, byrow = TRUE)
  df <- as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- as.character(seq_len(n_cols))
  rownames(df) <- LETTERS[seq_len(n_rows)]
  df
}

flatten_plate_df <- function(df) {
  as.vector(t(as.matrix(df)))
}

make_plate_factor_map <- function(fm_long, layout) {
  n_rows <- layout$n_rows[[1]]
  n_cols <- layout$n_cols[[1]]
  
  list(
    condition_df = make_plate_df(fm_long$condition, n_rows, n_cols, fill = ""),
    receptor_df  = make_plate_df(fm_long$receptor, n_rows, n_cols, fill = "none"),
    treatment_df = make_plate_df(fm_long$treatment, n_rows, n_cols, fill = "VEH")
  )
}

plate_to_long_factor_map <- function(condition_df, receptor_df, treatment_df) {
  conds <- flatten_plate_df(condition_df)
  recs  <- flatten_plate_df(receptor_df)
  trts  <- flatten_plate_df(treatment_df)
  
  keep <- trimws(conds) != ""
  
  tibble(
    condition = conds[keep],
    receptor = recs[keep],
    treatment = trts[keep]
  ) %>%
    mutate(
      receptor = purrr::map_chr(receptor, canonicalize_receptor_edit),
      treatment = if_else(trimws(treatment) == "", "VEH", trimws(treatment)),
      receptor = factor(receptor),
      treatment = factor(treatment)
    ) %>%
    arrange(condition)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Incucyte Multi-File Import + Channel Normalization"),
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
      h4("Factor editing"),
      actionButton("reset_factor_map", "Reset factor assignments"),
      actionButton("apply_factor_map", "Apply edited factor assignments", class = "btn-primary"),
      hr(),
      h4("Downloads"),
      downloadButton("download_factor_map", "Factor assignments (csv)"),
      downloadButton("download_norm", "Normalized by Passage (csv)"),
      downloadButton("download_prism", "Prism wide export (csv)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Plot", plotOutput("plot", height = 460)),
        tabPanel(
          "Factor editor",
          br(),
          uiOutput("plate_info"),
          tags$p("Edit receptor and treatment in plate layout. Condition header plate is read-only."),
          h4("Condition headers"),
          rHandsontableOutput("condition_plate", height = "240px"),
          br(),
          h4("Receptor"),
          rHandsontableOutput("receptor_plate", height = "240px"),
          br(),
          h4("Treatment"),
          rHandsontableOutput("treatment_plate", height = "240px")
        ),
        tabPanel("Join checks", verbatimTextOutput("join_checks")),
        tabPanel("File preview", tableOutput("preview_files")),
        tabPanel("Normalized", tableOutput("preview_norm")),
        tabPanel("Prism preview", tableOutput("preview_prism"))
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  
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
      file = input$files$name,
      path = input$files$datapath,
      channel = purrr::map_chr(seq_along(input$files$name), ~ get_chan(.x, input$files$name[.x]))
    )
  })
  
  output$preview_files <- renderTable({
    req(input$files)
    tibble(file = input$files$name, path = input$files$datapath) %>%
      mutate(preview = purrr::map(path, preview_file_lines, n = 10)) %>%
      select(file, preview) %>%
      tidyr::unnest(preview) %>%
      select(file, line_no, line)
  }, striped = TRUE)
  
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
          "Log2 ratio (log2(signal/control))" = "log2ratio",
          "OLS-adjusted (control-adjusted signal)" = "ols_adj"
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
  
  # ---------------------------
  # Imported data with auto-guessed factors
  # ---------------------------
  raw_long_auto <- eventReactive(input$run, {
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
    
    annot <- guess_receptor_treatment(unique(dat$condition))
    
    dat %>%
      left_join(annot, by = c("condition" = "condition_raw")) %>%
      mutate(
        receptor = forcats::fct_explicit_na(receptor, na_level = "none"),
        treatment = forcats::fct_explicit_na(treatment, na_level = "VEH")
      ) %>%
      relocate(receptor, treatment, .after = condition)
  }, ignoreInit = TRUE)
  
  factor_map_default <- reactive({
    req(raw_long_auto())
    raw_long_auto() %>%
      distinct(condition, receptor, treatment) %>%
      mutate(
        receptor = as.character(receptor),
        treatment = as.character(treatment)
      ) %>%
      arrange(condition)
  })
  
  plate_layout <- reactive({
    req(factor_map_default())
    detect_plate_layout(nrow(factor_map_default()))
  })
  
  output$plate_info <- renderUI({
    req(plate_layout())
    layout <- plate_layout()
    tags$p(
      tags$b("Detected plate layout: "),
      paste0(layout$plate_size[[1]], "-well equivalent (",
             layout$n_rows[[1]], " rows × ",
             layout$n_cols[[1]], " columns)")
    )
  })
  
  condition_plate_rv <- reactiveVal(NULL)
  receptor_plate_rv  <- reactiveVal(NULL)
  treatment_plate_rv <- reactiveVal(NULL)
  
  observeEvent(list(factor_map_default(), plate_layout()), {
    pm <- make_plate_factor_map(factor_map_default(), plate_layout())
    condition_plate_rv(pm$condition_df)
    receptor_plate_rv(pm$receptor_df)
    treatment_plate_rv(pm$treatment_df)
  })
  
  observeEvent(input$reset_factor_map, {
    req(factor_map_default(), plate_layout())
    pm <- make_plate_factor_map(factor_map_default(), plate_layout())
    condition_plate_rv(pm$condition_df)
    receptor_plate_rv(pm$receptor_df)
    treatment_plate_rv(pm$treatment_df)
  })
  
  output$condition_plate <- renderRHandsontable({
    req(condition_plate_rv())
    rhandsontable(
      condition_plate_rv(),
      rowHeaders = rownames(condition_plate_rv()),
      stretchH = "all",
      height = 220
    ) %>%
      hot_table(readOnly = TRUE, highlightCol = TRUE, highlightRow = TRUE)
  })
  
  output$receptor_plate <- renderRHandsontable({
    req(receptor_plate_rv())
    rhandsontable(
      receptor_plate_rv(),
      rowHeaders = rownames(receptor_plate_rv()),
      stretchH = "all",
      height = 220
    ) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })
  
  output$treatment_plate <- renderRHandsontable({
    req(treatment_plate_rv())
    rhandsontable(
      treatment_plate_rv(),
      rowHeaders = rownames(treatment_plate_rv()),
      stretchH = "all",
      height = 220
    ) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE)
  })
  
  observe({
    if (!is.null(input$receptor_plate)) {
      tbl <- hot_to_r(input$receptor_plate)
      if (!is.null(tbl)) {
        tbl <- as.data.frame(tbl, stringsAsFactors = FALSE, check.names = FALSE)
        rownames(tbl) <- rownames(receptor_plate_rv() %||% tbl)
        receptor_plate_rv(tbl)
      }
    }
  })
  
  observe({
    if (!is.null(input$treatment_plate)) {
      tbl <- hot_to_r(input$treatment_plate)
      if (!is.null(tbl)) {
        tbl <- as.data.frame(tbl, stringsAsFactors = FALSE, check.names = FALSE)
        rownames(tbl) <- rownames(treatment_plate_rv() %||% tbl)
        treatment_plate_rv(tbl)
      }
    }
  })
  
  edited_factor_map <- reactiveVal(NULL)
  
  observeEvent(input$apply_factor_map, {
    req(condition_plate_rv(), receptor_plate_rv(), treatment_plate_rv())
    
    df_long <- plate_to_long_factor_map(
      condition_df = condition_plate_rv(),
      receptor_df  = receptor_plate_rv(),
      treatment_df = treatment_plate_rv()
    )
    
    edited_factor_map(df_long)
  }, ignoreInit = TRUE)
  
  current_factor_map <- reactive({
    if (!is.null(edited_factor_map())) {
      edited_factor_map()
    } else {
      req(factor_map_default())
      factor_map_default() %>%
        mutate(
          receptor = factor(purrr::map_chr(receptor, canonicalize_receptor_edit)),
          treatment = factor(if_else(trimws(treatment) == "", "VEH", trimws(treatment)))
        )
    }
  })
  
  raw_long <- reactive({
    req(raw_long_auto(), current_factor_map())
    
    raw_long_auto() %>%
      select(-receptor, -treatment) %>%
      left_join(current_factor_map(), by = "condition") %>%
      mutate(
        receptor = forcats::fct_explicit_na(receptor, na_level = "none"),
        treatment = forcats::fct_explicit_na(treatment, na_level = "VEH")
      ) %>%
      relocate(receptor, treatment, .after = condition)
  })
  
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
    
    factor_summary <- rl %>%
      distinct(condition, receptor, treatment) %>%
      summarize(
        n_conditions = n(),
        receptor_levels = paste(sort(unique(as.character(receptor))), collapse = ", "),
        treatment_levels = paste(sort(unique(as.character(treatment))), collapse = ", ")
      )
    
    list(
      files_loaded = files_loaded,
      per_channel_summary = per_channel,
      factor_assignment_summary = factor_summary,
      note = "Factor editor starts from auto-guessed values, auto-detects plate size, and displays conditions in plate-equivalent layout."
    )
  })
  
  wide_joined_passage <- reactive({
    req(raw_long())
    raw_long() %>%
      select(passage, channel, elapsed, condition, receptor, treatment, value) %>%
      group_by(passage, channel, elapsed, condition, receptor, treatment) %>%
      summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
      pivot_wider(names_from = channel, values_from = value)
  })
  
  normalized_passage <- reactive({
    req(wide_joined_passage())
    w <- wide_joined_passage()
    
    method <- input$norm_method %||% "ratio"
    sig <- input$signal_channel
    ctl <- input$control_channel
    
    if (is.null(sig) || is.null(ctl) || !(sig %in% names(w)) || !(ctl %in% names(w))) {
      return(w %>% mutate(value_norm = NA_real_))
    }
    
    out <- w
    
    if (method %in% c("none", "ratio", "log2ratio")) {
      out <- out %>%
        mutate(
          value_norm = dplyr::case_when(
            method == "none" ~ .data[[sig]],
            method == "ratio" ~ .data[[sig]] / .data[[ctl]],
            method == "log2ratio" ~ log2(.data[[sig]] / .data[[ctl]]),
            TRUE ~ NA_real_
          )
        )
    } else if (method == "ols_adj") {
      out <- out %>%
        group_by(passage, condition) %>%
        group_modify(~ ols_control_adjust(.x, sig_col = sig, ctl_col = ctl)) %>%
        ungroup()
    } else {
      out <- out %>% mutate(value_norm = NA_real_)
    }
    
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
  }, striped = TRUE)
  
  stats_long <- reactive({
    req(normalized_passage())
    normalized_passage() %>%
      mutate(
        receptor = forcats::fct_explicit_na(receptor, na_level = "none"),
        treatment = forcats::fct_explicit_na(treatment, na_level = "VEH")
      ) %>%
      select(passage, elapsed, condition, receptor, treatment, value_norm) %>%
      arrange(receptor, treatment, passage, condition, elapsed)
  })
  
  prism_wide <- reactive({
    req(stats_long())
    df <- stats_long() %>%
      group_by(elapsed, receptor, treatment, passage) %>%
      summarize(value = mean(value_norm, na.rm = TRUE), .groups = "drop") %>%
      mutate(col_key = paste(receptor, treatment, passage, sep = "__"))
    
    wide <- df %>%
      select(elapsed, col_key, value) %>%
      pivot_wider(names_from = col_key, values_from = value)
    
    col_meta <- df %>%
      distinct(col_key, receptor, treatment, passage) %>%
      arrange(receptor, treatment, passage)
    
    wide %>% select(elapsed, any_of(col_meta$col_key))
  })
  
  output$preview_prism <- renderTable({
    req(prism_wide())
    head(prism_wide(), 20)
  }, striped = TRUE)
  
  output$plot <- renderPlot({
    req(stats_long())
    df <- stats_long()
    
    if (all(is.na(df$value_norm))) {
      plot.new()
      text(0.5, 0.5, "No normalized values to plot yet.\nUpload files, assign channels, choose signal/control, click 'Import + Process'.")
      return()
    }
    
    ggplot(df, aes(x = elapsed, y = value_norm, group = condition, color = receptor)) +
      geom_line(alpha = 0.6, na.rm = TRUE) +
      facet_wrap(~ passage) +
      labs(
        x = "Elapsed time",
        y = "Value (normalized)",
        title = "Normalized trajectories (facet by Passage; color by receptor)"
      ) +
      theme_minimal()
  })
  
  output$download_factor_map <- downloadHandler(
    filename = function() paste0("incucyte_factor_assignments_", Sys.Date(), ".csv"),
    content = function(file) {
      req(current_factor_map())
      write_csv(current_factor_map(), file)
    }
  )
  
  output$download_norm <- downloadHandler(
    filename = function() paste0("incucyte_normalized_by_passage_", Sys.Date(), ".csv"),
    content = function(file) {
      req(normalized_passage())
      write_csv(normalized_passage(), file)
    }
  )
  
  output$download_prism <- downloadHandler(
    filename = function() paste0("incucyte_prism_wide_1row_per_time_", Sys.Date(), ".csv"),
    content = function(file) {
      req(prism_wide())
      write_csv(prism_wide(), file)
    }
  )
}

shinyApp(ui, server)