# app.R
# Incucyte Multi-File Import + Channel Normalization
# - Auto-parses receptor + treatment from condition headers
# - Extracts well IDs like A3 from headers and shows them in the editor
# - Matching across files/channels is based on parsed factors, not raw condition strings
# - Factor editor uses rhandsontable in LONG format with sortable columns
# - Passage is editable in the same factor editor table, so mixed-passage plates are supported
# - Receptor canonical names after parsing:
#     era/er -> ER_a
#     erb    -> ER_b
#     pgr/pr -> PR
#     pra    -> PR_a
#     prb    -> PR_b
#     ar     -> AR
#     gr     -> GR
#     mr     -> MR
# - Receptor combinations are canonicalized, so equivalent combos always sort identically
# - Treatment combinations are canonicalized, so equivalent combos always sort identically
# - Factor editor includes:
#     * file
#     * well_id
#     * condition
#     * n_reps
#     * passage
#     * receptor
#     * treatment
# - Plot tab includes:
#     * filtered timecourse plot
#     * combined preview AUC plot for the current filter window
# - Prism export uses AUC values over the currently selected plot/filter window
#   formatted as rows = receptor, columns grouped by treatment, one column per replicate

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

extract_well_id <- function(x) {
  m <- str_match(x, "\\(([A-H][0-9]{1,2})\\)")
  out <- m[, 2]
  out[is.na(out)] <- ""
  out
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
      elapsed  = suppressWarnings(as.numeric(elapsed))
    ) %>%
    pivot_longer(
      cols = -c(datetime, elapsed),
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(
      condition = as.character(condition),
      value     = suppressWarnings(as.numeric(value)),
      well_id   = extract_well_id(condition)
    )
}

# ---------------------------
# Header cleaning / parsing
# ---------------------------
clean_condition_header <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("\\([a-h][0-9]{1,2}\\)", " ") %>%           # remove well IDs like (A1)
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*ul\\b", " ") %>%    # remove dispense volumes
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*mg/ml\\b", " ") %>% # remove stock concentration notes
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

canonicalize_receptor_combo <- function(x) {
  if (is.na(x) || x == "" || tolower(x) == "none") return("none")
  
  parts <- str_split(as.character(x), "\\s*\\+\\s*", simplify = FALSE)[[1]]
  parts <- str_squish(parts)
  parts <- parts[parts != ""]
  parts <- unique(parts)
  
  canon_part <- function(p) {
    z <- p %>%
      str_to_lower() %>%
      str_replace_all("[α]", "a") %>%
      str_replace_all("[β]", "b") %>%
      str_replace_all("[^a-z0-9_ ]", " ") %>%
      str_squish()
    
    dplyr::case_when(
      z %in% c("era", "er", "er a", "esr1", "er1", "er_a") ~ "ER_a",
      z %in% c("erb", "er b", "esr2", "er2", "er_b") ~ "ER_b",
      z %in% c("pgr", "pr") ~ "PR",
      z %in% c("pra", "pr a", "pr_a") ~ "PR_a",
      z %in% c("prb", "pr b", "pr_b") ~ "PR_b",
      z %in% c("ar") ~ "AR",
      z %in% c("gr") ~ "GR",
      z %in% c("mr") ~ "MR",
      z %in% c("", "none", "na") ~ "none",
      TRUE ~ str_trim(p)
    )
  }
  
  parts <- purrr::map_chr(parts, canon_part)
  parts <- unique(parts[parts != "none" & parts != ""])
  
  # if subtype present, drop generic PR
  if (any(c("PR_a", "PR_b") %in% parts)) parts <- setdiff(parts, "PR")
  
  ord <- c("ER_a", "ER_b", "PR", "PR_a", "PR_b", "AR", "GR", "MR")
  parts <- parts[order(match(parts, ord, nomatch = 999), parts)]
  
  if (length(parts) == 0) "none" else paste(parts, collapse = " + ")
}

canonicalize_treatment_combo <- function(x) {
  z <- str_trim(as.character(x))
  if (is.na(z) || z == "" || toupper(z) == "VEH") return("VEH")
  
  parts <- str_split(z, "\\s*\\+\\s*", simplify = FALSE)[[1]]
  parts <- toupper(str_squish(parts))
  parts <- parts[parts != ""]
  parts <- unique(parts)
  
  get_ligand <- function(s) {
    m <- str_match(s, "\\b(E2|P4|DHT)\\b")
    lig <- m[, 2]
    ifelse(is.na(lig), "ZZZ", lig)
  }
  
  get_dose <- function(s) {
    m <- str_match(s, "^([0-9]+\\.?[0-9]*)")
    dose <- suppressWarnings(as.numeric(m[, 2]))
    ifelse(is.na(dose), Inf, dose)
  }
  
  ligand_order <- c("E2", "P4", "DHT", "ZZZ")
  ord_lig <- match(get_ligand(parts), ligand_order)
  ord_dose <- get_dose(parts)
  
  parts <- parts[order(ord_lig, ord_dose, parts)]
  paste(parts, collapse = " + ")
}

extract_receptors <- function(x) {
  s <- clean_condition_header(x)
  recs <- c()
  
  if (str_detect(s, "\\berb\\b|\\ber b\\b|\\besr2\\b|\\ber2\\b")) recs <- c(recs, "ER_b")
  if (str_detect(s, "\\ber\\b|\\bera\\b|\\ber a\\b|\\besr1\\b|\\ber1\\b")) recs <- c(recs, "ER_a")
  if (str_detect(s, "\\bpra\\b|\\bpr a\\b")) recs <- c(recs, "PR_a")
  if (str_detect(s, "\\bprb\\b|\\bpr b\\b")) recs <- c(recs, "PR_b")
  if (str_detect(s, "\\bpgr\\b|\\bpr\\b")) recs <- c(recs, "PR")
  if (str_detect(s, "\\bar\\b|\\bandrogen receptor\\b|\\bnr3c4\\b")) recs <- c(recs, "AR")
  if (str_detect(s, "\\bgr\\b|\\bglucocorticoid receptor\\b|\\bnr3c1\\b")) recs <- c(recs, "GR")
  if (str_detect(s, "\\bmr\\b|\\bmineralocorticoid receptor\\b|\\bnr3c2\\b")) recs <- c(recs, "MR")
  
  canonicalize_receptor_combo(paste(unique(recs), collapse = " + "))
}

extract_treatment <- function(x) {
  s <- clean_condition_header(x)
  veh_only <- str_detect(s, "\\bveh\\b|\\bvehicle\\b|\\be2oh\\b|\\bethanol\\b|\\bdmso\\b")
  
  tokens <- str_split(s, "\\s+", simplify = TRUE)
  tokens <- tokens[tokens != ""]
  
  is_num <- function(z) str_detect(z, "^[0-9]+\\.?[0-9]*$")
  is_unit <- function(z) str_detect(z, regex("^(pm|nm|um|µm|mm)$", ignore_case = TRUE))
  is_ligand <- function(z) str_detect(z, regex("^(e2|p4|dht)$", ignore_case = TRUE))
  
  normalize_unit <- function(z) {
    z <- toupper(z)
    if (z == "µM") "uM" else z
  }
  
  normalize_ligand <- function(z) toupper(z)
  
  hits <- character(0)
  used <- rep(FALSE, length(tokens))
  
  i <- 1
  while (i <= length(tokens)) {
    if (used[i]) {
      i <- i + 1
      next
    }
    
    # ligand-first: E2 30 nm
    if (i + 2 <= length(tokens)) {
      if (!used[i] && !used[i + 1] && !used[i + 2] &&
          is_ligand(tokens[i]) &&
          is_num(tokens[i + 1]) &&
          is_unit(tokens[i + 2])) {
        
        hits <- c(
          hits,
          paste0(
            tokens[i + 1],
            normalize_unit(tokens[i + 2]),
            " ",
            normalize_ligand(tokens[i])
          )
        )
        used[i:(i + 2)] <- TRUE
        i <- i + 3
        next
      }
    }
    
    # amount-first: 30 nm E2
    if (i + 2 <= length(tokens)) {
      if (!used[i] && !used[i + 1] && !used[i + 2] &&
          is_num(tokens[i]) &&
          is_unit(tokens[i + 1]) &&
          is_ligand(tokens[i + 2])) {
        
        hits <- c(
          hits,
          paste0(
            tokens[i],
            normalize_unit(tokens[i + 1]),
            " ",
            normalize_ligand(tokens[i + 2])
          )
        )
        used[i:(i + 2)] <- TRUE
        i <- i + 3
        next
      }
    }
    
    i <- i + 1
  }
  
  hits <- unique(hits)
  
  if (length(hits) > 0) return(canonicalize_treatment_combo(paste(hits, collapse = " + ")))
  if (veh_only) return("VEH")
  "VEH"
}

canonicalize_receptor_edit <- function(x) {
  canonicalize_receptor_combo(x)
}

canonicalize_treatment_edit <- function(x) {
  z <- str_trim(as.character(x))
  if (is.na(z) || z == "") return("VEH")
  canonicalize_treatment_combo(toupper(z))
}

canonicalize_passage_edit <- function(x) {
  z <- str_trim(as.character(x))
  if (is.na(z) || z == "") "Passage_NA" else z
}

normalize_factor_key <- function(x) {
  x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("\\s+", " ") %>%
    str_squish()
}

make_factor_key <- function(receptor, treatment) {
  paste(normalize_factor_key(receptor), normalize_factor_key(treatment), sep = " || ")
}

guess_receptor_treatment <- function(condition_vec) {
  tibble(condition_raw = condition_vec) %>%
    mutate(
      receptor = purrr::map_chr(condition_raw, extract_receptors),
      treatment = purrr::map_chr(condition_raw, extract_treatment),
      receptor = factor(receptor),
      treatment = factor(treatment)
    )
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

auc_trapz <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 2) return(NA_real_)
  
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
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
      h4("Editing"),
      actionButton("reset_editor", "Reset editor"),
      actionButton("apply_editor", "Apply edits", class = "btn-primary"),
      hr(),
      h4("Downloads"),
      downloadButton("download_editor", "Editor table (csv)"),
      downloadButton("download_norm", "Normalized by Passage (csv)"),
      downloadButton("download_prism", "Prism AUC export (csv)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Plot",
          plotOutput("plot", height = 420),
          br(),
          plotOutput("auc_plot", height = 340),
          br(),
          fluidRow(
            column(4, uiOutput("plot_time_ui")),
            column(4, uiOutput("plot_receptor_ui")),
            column(4, uiOutput("plot_treatment_ui"))
          )
        ),
        tabPanel(
          "Factor editor",
          br(),
          tags$p("Edit passage, receptor, and treatment. Columns are sortable. The table supports copy/paste and drag-fill. Click 'Apply edits' when done."),
          rHandsontableOutput("editor_table", height = "560px")
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
        "Also baseline-normalize within Passage+Factor combination (divide by first non-NA timepoint)",
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
          passage = passage_lab,
          vessel_name = meta$vessel_name[[1]] %||% NA_character_,
          metric      = meta$metric[[1]] %||% NA_character_,
          cell_type   = meta$cell_type[[1]] %||% NA_character_,
          analysis    = meta$analysis[[1]] %||% NA_character_,
          condition_id = paste(file, condition, sep = " || ")
        )
    }) %>%
      select(condition_id, file, well_id, channel, passage, vessel_name, metric, cell_type, analysis,
             datetime, elapsed, condition, value)
    
    annot <- guess_receptor_treatment(unique(dat$condition))
    
    dat %>%
      left_join(annot, by = c("condition" = "condition_raw")) %>%
      mutate(
        receptor = forcats::fct_explicit_na(receptor, na_level = "none"),
        treatment = forcats::fct_explicit_na(treatment, na_level = "VEH"),
        receptor = factor(purrr::map_chr(as.character(receptor), canonicalize_receptor_combo)),
        treatment = factor(purrr::map_chr(as.character(treatment), canonicalize_treatment_combo)),
        factor_key = make_factor_key(receptor, treatment)
      ) %>%
      relocate(receptor, treatment, factor_key, .after = condition)
  }, ignoreInit = TRUE)
  
  editor_default <- reactive({
    req(raw_long_auto())
    
    rep_counts <- raw_long_auto() %>%
      distinct(condition, file) %>%
      count(condition, name = "n_reps")
    
    raw_long_auto() %>%
      distinct(condition_id, file, well_id, condition, passage, receptor, treatment, factor_key) %>%
      left_join(rep_counts, by = "condition") %>%
      mutate(
        well_id    = as.character(well_id),
        passage    = as.character(passage),
        receptor   = as.character(receptor),
        treatment  = as.character(treatment),
        factor_key = as.character(factor_key),
        n_reps     = as.integer(n_reps)
      ) %>%
      arrange(file, well_id, condition)
  })
  
  editor_rv <- reactiveVal(NULL)
  
  observeEvent(editor_default(), {
    editor_rv(editor_default())
  })
  
  observeEvent(input$reset_editor, {
    req(editor_default())
    editor_rv(editor_default())
  })
  
  output$editor_table <- renderRHandsontable({
    req(editor_rv())
    df <- editor_rv() %>%
      select(file, well_id, condition, n_reps, passage, receptor, treatment)
    
    rhandsontable(
      df,
      rowHeaders = NULL,
      stretchH = "all",
      height = 540
    ) %>%
      hot_col("file", readOnly = TRUE) %>%
      hot_col("well_id", readOnly = TRUE) %>%
      hot_col("condition", readOnly = TRUE) %>%
      hot_col("n_reps", readOnly = TRUE) %>%
      hot_table(
        highlightCol = TRUE,
        highlightRow = TRUE,
        columnSorting = TRUE,
        manualColumnMove = TRUE
      )
  })
  
  observe({
    if (!is.null(input$editor_table)) {
      tbl <- hot_to_r(input$editor_table)
      if (!is.null(tbl) && !is.null(editor_rv())) {
        tbl <- as_tibble(tbl) %>%
          mutate(
            file = as.character(file),
            well_id = as.character(well_id),
            condition = as.character(condition),
            n_reps = suppressWarnings(as.integer(n_reps)),
            passage = as.character(passage),
            receptor = as.character(receptor),
            treatment = as.character(treatment)
          )
        
        key_map <- editor_rv() %>%
          transmute(condition_id, file, well_id, condition) %>%
          distinct()
        
        updated <- tbl %>%
          left_join(key_map, by = c("file", "well_id", "condition")) %>%
          select(condition_id, everything())
        
        editor_rv(updated)
      }
    }
  })
  
  edited_map <- reactiveVal(NULL)
  
  observeEvent(input$apply_editor, {
    req(editor_rv())
    
    df_edit <- editor_rv() %>%
      mutate(
        well_id = as.character(well_id),
        passage = purrr::map_chr(passage, canonicalize_passage_edit),
        receptor = purrr::map_chr(receptor, canonicalize_receptor_edit),
        treatment = purrr::map_chr(treatment, canonicalize_treatment_edit),
        factor_key = make_factor_key(receptor, treatment),
        receptor = factor(receptor),
        treatment = factor(treatment),
        passage = factor(passage)
      ) %>%
      arrange(file, well_id, condition)
    
    edited_map(df_edit)
  }, ignoreInit = TRUE)
  
  current_editor_map <- reactive({
    if (!is.null(edited_map())) {
      edited_map()
    } else {
      req(editor_default())
      editor_default() %>%
        mutate(
          well_id = as.character(well_id),
          passage = factor(purrr::map_chr(passage, canonicalize_passage_edit)),
          receptor = factor(purrr::map_chr(receptor, canonicalize_receptor_edit)),
          treatment = factor(purrr::map_chr(treatment, canonicalize_treatment_edit)),
          factor_key = make_factor_key(receptor, treatment)
        )
    }
  })
  
  raw_long <- reactive({
    req(raw_long_auto(), current_editor_map())
    
    raw_long_auto() %>%
      select(-passage, -receptor, -treatment, -factor_key) %>%
      left_join(
        current_editor_map() %>%
          select(condition_id, passage, receptor, treatment, factor_key),
        by = "condition_id"
      ) %>%
      mutate(
        passage = factor(as.character(passage)),
        receptor = forcats::fct_explicit_na(receptor, na_level = "none"),
        treatment = forcats::fct_explicit_na(treatment, na_level = "VEH"),
        receptor = factor(purrr::map_chr(as.character(receptor), canonicalize_receptor_combo)),
        treatment = factor(purrr::map_chr(as.character(treatment), canonicalize_treatment_combo)),
        factor_key = make_factor_key(receptor, treatment)
      ) %>%
      relocate(well_id, passage, receptor, treatment, factor_key, .after = condition)
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
        n_factor_keys = n_distinct(factor_key),
        n_times = n_distinct(elapsed),
        .groups = "drop"
      )
    
    factor_summary <- rl %>%
      distinct(condition_id, condition, well_id, passage, receptor, treatment, factor_key) %>%
      summarize(
        n_rows = n(),
        n_unique_factor_keys = n_distinct(factor_key),
        receptor_levels = paste(sort(unique(as.character(receptor))), collapse = ", "),
        treatment_levels = paste(sort(unique(as.character(treatment))), collapse = ", "),
        passage_levels = paste(sort(unique(as.character(passage))), collapse = ", ")
      )
    
    list(
      files_loaded = files_loaded,
      per_channel_summary = per_channel,
      factor_assignment_summary = factor_summary,
      note = "Files/channels are matched by parsed factor combination (receptor + treatment). Passage is editable per imported condition row, mixed-passage plates are supported, and well IDs are extracted from headers when present."
    )
  })
  
  wide_joined_passage <- reactive({
    req(raw_long())
    raw_long() %>%
      select(passage, channel, elapsed, receptor, treatment, factor_key, value) %>%
      group_by(passage, channel, elapsed, receptor, treatment, factor_key) %>%
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
        group_by(passage, factor_key) %>%
        group_modify(~ ols_control_adjust(.x, sig_col = sig, ctl_col = ctl)) %>%
        ungroup()
    } else {
      out <- out %>% mutate(value_norm = NA_real_)
    }
    
    if (isTRUE(input$baseline_norm)) {
      out <- out %>%
        group_by(passage, factor_key) %>%
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
      select(passage, elapsed, receptor, treatment, factor_key, value_norm) %>%
      arrange(receptor, treatment, passage, elapsed)
  })
  
  output$plot_time_ui <- renderUI({
    req(stats_long())
    df <- stats_long() %>% filter(is.finite(elapsed))
    
    if (nrow(df) == 0) return(NULL)
    
    rng <- range(df$elapsed, na.rm = TRUE)
    step_val <- max((rng[2] - rng[1]) / 100, .Machine$double.eps)
    
    sliderInput(
      "plot_time_range",
      "Elapsed time range",
      min = floor(rng[1]),
      max = ceiling(rng[2]),
      value = c(floor(rng[1]), ceiling(rng[2])),
      step = step_val
    )
  })
  
  output$plot_receptor_ui <- renderUI({
    req(stats_long())
    levs <- stats_long() %>%
      pull(receptor) %>%
      as.character() %>%
      unique() %>%
      sort()
    
    checkboxGroupInput(
      "plot_receptors",
      "Receptors to show",
      choices = levs,
      selected = levs
    )
  })
  
  output$plot_treatment_ui <- renderUI({
    req(stats_long())
    levs <- stats_long() %>%
      pull(treatment) %>%
      as.character() %>%
      unique() %>%
      sort()
    
    checkboxGroupInput(
      "plot_treatments",
      "Treatments to show",
      choices = levs,
      selected = levs
    )
  })
  
  filtered_stats_long <- reactive({
    req(stats_long())
    df <- stats_long()
    
    if (!is.null(input$plot_time_range)) {
      df <- df %>%
        filter(
          elapsed >= input$plot_time_range[1],
          elapsed <= input$plot_time_range[2]
        )
    }
    
    if (!is.null(input$plot_receptors) && length(input$plot_receptors) > 0) {
      df <- df %>%
        filter(as.character(receptor) %in% input$plot_receptors)
    }
    
    if (!is.null(input$plot_treatments) && length(input$plot_treatments) > 0) {
      df <- df %>%
        filter(as.character(treatment) %in% input$plot_treatments)
    }
    
    df
  })
  
  auc_preview <- reactive({
    req(filtered_stats_long())
    filtered_stats_long() %>%
      group_by(receptor, treatment, passage) %>%
      summarize(
        auc = auc_trapz(elapsed, value_norm),
        .groups = "drop"
      )
  })
  
  prism_wide <- reactive({
    req(auc_preview())
    auc_df <- auc_preview() %>%
      arrange(receptor, treatment, passage)
    
    auc_df <- auc_df %>%
      group_by(receptor, treatment) %>%
      mutate(rep_idx = row_number()) %>%
      ungroup()
    
    col_template <- auc_df %>%
      count(treatment, name = "n_rep") %>%
      group_by(treatment) %>%
      summarise(max_rep = max(n_rep), .groups = "drop") %>%
      mutate(col_keys = purrr::map2(treatment, max_rep, ~ paste0(.x, "__rep", seq_len(.y)))) %>%
      pull(col_keys) %>%
      unlist()
    
    wide <- auc_df %>%
      mutate(col_key = paste0(treatment, "__rep", rep_idx)) %>%
      select(receptor, col_key, auc) %>%
      pivot_wider(names_from = col_key, values_from = auc)
    
    wide %>%
      select(receptor, any_of(col_template))
  })
  
  output$preview_prism <- renderTable({
    req(prism_wide())
    prism_wide()
  }, striped = TRUE)
  
  output$plot <- renderPlot({
    req(filtered_stats_long())
    df <- filtered_stats_long()
    
    if (all(is.na(df$value_norm)) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data remain after applying the current plot filters.")
      return()
    }
    
    ggplot(df, aes(x = elapsed, y = value_norm, group = factor_key, color = receptor)) +
      geom_line(alpha = 0.6, na.rm = TRUE) +
      facet_wrap(~ passage) +
      labs(
        x = "Elapsed time",
        y = "Value (normalized)",
        title = "Normalized trajectories (facet by Passage; color by receptor)"
      ) +
      theme_minimal()
  })
  
  output$auc_plot <- renderPlot({
    req(auc_preview())
    df <- auc_preview()
    
    if (nrow(df) == 0 || all(!is.finite(df$auc))) {
      plot.new()
      text(0.5, 0.5, "No AUC values available for the current filters.")
      return()
    }
    
    dodge <- position_dodge(width = 0.6)
    
    ggplot(df, aes(x = treatment, y = auc, color = receptor)) +
      geom_point(
        position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.6),
        alpha = 0.75,
        size = 2
      ) +
      stat_summary(aes(group = receptor), fun = mean, geom = "line", position = dodge, linewidth = 0.8) +
      stat_summary(aes(group = receptor), fun = mean, geom = "point", position = dodge, size = 2.5) +
      labs(
        x = "Treatment",
        y = "AUC",
        title = "Combined AUC preview for current filter window"
      ) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$download_editor <- downloadHandler(
    filename = function() paste0("incucyte_editor_table_", Sys.Date(), ".csv"),
    content = function(file) {
      req(current_editor_map())
      write_csv(current_editor_map(), file)
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
    filename = function() paste0("incucyte_prism_auc_", Sys.Date(), ".csv"),
    content = function(file) {
      req(prism_wide())
      
      wide <- prism_wide()
      value_cols <- names(wide)[-1]
      
      treatment_header <- c(
        "",
        stringr::str_replace(value_cols, "__rep\\d+$", "")
      )
      
      rep_header <- c(
        "receptor",
        stringr::str_extract(value_cols, "rep\\d+$")
      )
      
      out <- rbind(
        treatment_header,
        rep_header,
        as.matrix(wide)
      )
      
      write.table(
        out,
        file = file,
        sep = ",",
        row.names = FALSE,
        col.names = FALSE,
        quote = TRUE,
        na = ""
      )
    }
  )
}

shinyApp(ui, server)