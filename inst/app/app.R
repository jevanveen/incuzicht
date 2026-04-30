# app.R
# Incucyte Multi-File Import + Channel Normalization

library(shiny)
library(tidyverse)
library(readr)
library(stringr)
library(rhandsontable)
library(broom)
library(DT)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------
# Low-level helpers
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

read_incucyte_raw_table <- function(path) {
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
  dat
}

# ---------------------------
# Canonicalization
# ---------------------------
canonicalize_receptor_combo <- function(x) {
  if (length(x) != 1) return("none")
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
  
  if (any(c("PR_a", "PR_b") %in% parts)) parts <- setdiff(parts, "PR")
  
  ord <- c("ER_a", "ER_b", "PR", "PR_a", "PR_b", "AR", "GR", "MR")
  parts <- parts[order(match(parts, ord, nomatch = 999), parts)]
  
  if (length(parts) == 0) "none" else paste(parts, collapse = " + ")
}

canonicalize_treatment_combo <- function(x) {
  z <- str_trim(as.character(x))
  if (length(z) != 1 || is.na(z) || z == "" || toupper(z) == "VEH") return("VEH")
  
  parts <- str_split(z, "\\s*\\+\\s*", simplify = FALSE)[[1]]
  parts <- toupper(str_squish(parts))
  parts <- parts[parts != ""]
  parts <- unique(parts)
  
  get_ligand <- function(s) {
    m <- str_match(s, "\\b(E2|P4|DHT|4-OHT|DEX|CORT|RU-486)\\b")
    lig <- m[, 2]
    ifelse(is.na(lig), "ZZZ", lig)
  }
  
  get_dose <- function(s) {
    m <- str_match(s, "^([0-9]+\\.?[0-9]*)")
    dose <- suppressWarnings(as.numeric(m[, 2]))
    ifelse(is.na(dose), Inf, dose)
  }
  
  ligand_order <- c("E2", "P4", "DHT", "4-OHT", "RU-486", "DEX", "CORT", "ZZZ")
  ord_lig <- match(get_ligand(parts), ligand_order)
  ord_dose <- get_dose(parts)
  
  parts <- parts[order(ord_lig, ord_dose, parts)]
  paste(parts, collapse = " + ")
}

canonicalize_receptor_edit <- function(x) canonicalize_receptor_combo(x)

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

# ---------------------------
# Parsing
# ---------------------------
clean_condition_header <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("\\([a-h][0-9]{1,2}\\)", " ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*ul\\b", " ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*mg/ml\\b", " ") %>%
    str_replace_all("([0-9]+\\.?[0-9]*)(pm|nm|um|µm|mm)\\b", "\\1 \\2") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s*\\+\\s*", " + ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

extract_receptors <- function(x) {
  s <- clean_condition_header(x)
  recs <- c()
  
  if (str_detect(s, "\\berb\\b|\\ber b\\b|\\besr2\\b|\\ber2\\b")) recs <- c(recs, "ER_b")
  if (str_detect(s, "\\ber\\b|\\bera\\b|\\ber a\\b|\\besr1\\b|\\ber1\\b|\\bnoer\\b")) {
    if (str_detect(s, "\\bnoer\\b")) recs <- c(recs, "NOER") else recs <- c(recs, "ER_a")
  }
  if (str_detect(s, "\\bpra\\b|\\bpr a\\b")) recs <- c(recs, "PR_a")
  if (str_detect(s, "\\bprb\\b|\\bpr b\\b")) recs <- c(recs, "PR_b")
  if (str_detect(s, "\\bpgr\\b|\\bpr\\b")) recs <- c(recs, "PR")
  if (str_detect(s, "\\bar\\b|\\bandrogen receptor\\b|\\bnr3c4\\b")) recs <- c(recs, "AR")
  if (str_detect(s, "\\bgr\\b|\\bglucocorticoid receptor\\b|\\bnr3c1\\b")) recs <- c(recs, "GR")
  if (str_detect(s, "\\bmr\\b|\\bmineralocorticoid receptor\\b|\\bnr3c2\\b")) recs <- c(recs, "MR")
  
  if (length(recs) == 0) return("none")
  
  canonical <- recs[recs %in% c("ER_a", "ER_b", "PR", "PR_a", "PR_b", "AR", "GR", "MR")]
  custom <- setdiff(unique(recs), canonical)
  
  canonical <- canonicalize_receptor_combo(paste(canonical, collapse = " + "))
  canon_parts <- if (canonical == "none") character(0) else str_split(canonical, "\\s*\\+\\s*", simplify = FALSE)[[1]]
  
  final_parts <- unique(c(canon_parts, custom))
  if (length(final_parts) == 0) "none" else paste(final_parts, collapse = " + ")
}

extract_treatment <- function(x) {
  s <- clean_condition_header(x)
  veh_only <- str_detect(s, "\\bveh\\b|\\bvehicle\\b|\\be2oh\\b|\\bethanol\\b|\\bdmso\\b")
  
  tokens <- str_split(s, "\\s+", simplify = TRUE)
  tokens <- tokens[tokens != ""]
  
  is_num <- function(z) str_detect(z, "^[0-9]+\\.?[0-9]*$")
  is_unit <- function(z) str_detect(z, regex("^(pm|nm|um|µm|mm)$", ignore_case = TRUE))
  is_ligand <- function(z) str_detect(z, regex("^(e2|p4|dht|4-oht|dex|cort|ru-486)$", ignore_case = TRUE))
  
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
    
    if (i + 2 <= length(tokens)) {
      if (!used[i] && !used[i + 1] && !used[i + 2] &&
          is_ligand(tokens[i]) &&
          is_num(tokens[i + 1]) &&
          is_unit(tokens[i + 2])) {
        hits <- c(hits, paste0(tokens[i + 1], normalize_unit(tokens[i + 2]), " ", normalize_ligand(tokens[i])))
        used[i:(i + 2)] <- TRUE
        i <- i + 3
        next
      }
    }
    
    if (i + 2 <= length(tokens)) {
      if (!used[i] && !used[i + 1] && !used[i + 2] &&
          is_num(tokens[i]) &&
          is_unit(tokens[i + 1]) &&
          is_ligand(tokens[i + 2])) {
        hits <- c(hits, paste0(tokens[i], normalize_unit(tokens[i + 1]), " ", normalize_ligand(tokens[i + 2])))
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
  
  s2 <- s %>%
    str_replace_all("\\b(era|er a|erb|er b|er|esr1|esr2|pra|pr a|prb|pr b|pgr|pr|ar|gr|mr|androgen receptor|glucocorticoid receptor|mineralocorticoid receptor)\\b", " ") %>%
    str_replace_all("\\s*\\+\\s*", " + ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
  
  s2 <- str_replace_all(s2, "\\bnone\\b", "")
  s2 <- str_replace_all(s2, "\\s+", " ")
  s2 <- str_trim(s2)
  
  if (s2 == "" || s2 == "+") return("VEH")
  toupper(s2)
}

parse_comma_header <- function(header) {
  if (str_detect(header, "Std Err")) {
    return(tibble(
      condition = header,
      well_id = "",
      receptor = NA_character_,
      treatment = NA_character_,
      replicate_id = NA_character_,
      drop_me = TRUE
    ))
  }
  
  parts <- str_split(header, ",")[[1]]
  treatment_part <- if (length(parts) >= 1) parts[1] else header
  rhs <- if (length(parts) >= 2) paste(parts[-1], collapse = ",") else ""
  
  rep_id <- str_match(rhs, "\\((\\d+)\\)")[, 2]
  receptor_token <- str_match(rhs, "\\+\\s*([A-Za-z0-9_\\-]+)\\s*\\(")[, 2]
  
  receptor <- ifelse(
    is.na(receptor_token),
    extract_receptors(rhs),
    canonicalize_receptor_edit(receptor_token)
  )
  
  tibble(
    condition = header,
    well_id = "",
    receptor = receptor,
    treatment = extract_treatment(treatment_part),
    replicate_id = ifelse(is.na(rep_id), "", rep_id),
    drop_me = FALSE
  )
}

guess_receptor_treatment <- function(condition_vec) {
  tibble(condition_raw = condition_vec) %>%
    mutate(
      receptor = purrr::map_chr(condition_raw, extract_receptors),
      treatment = purrr::map_chr(condition_raw, extract_treatment)
    )
}

read_incucyte_long <- function(path, drop_stderr = TRUE) {
  dat <- read_incucyte_raw_table(path)
  if (drop_stderr) dat <- dat %>% select(-matches("Std Err"))
  
  cond_names <- names(dat)[-(1:2)]
  is_comma_schema <- any(str_detect(cond_names, ","))
  
  long <- dat %>%
    pivot_longer(
      cols = -c(datetime, elapsed),
      names_to = "condition",
      values_to = "value"
    ) %>%
    mutate(
      condition = as.character(condition),
      value = suppressWarnings(as.numeric(value)),
      elapsed = suppressWarnings(as.numeric(elapsed)),
      datetime = as.character(datetime)
    )
  
  if (is_comma_schema) {
    parsed <- purrr::map_dfr(unique(long$condition), parse_comma_header)
    
    long %>%
      left_join(parsed, by = "condition") %>%
      filter(!drop_me) %>%
      mutate(
        well_id = "",
        receptor = if_else(is.na(receptor), "none", receptor),
        treatment = if_else(is.na(treatment), "VEH", treatment)
      ) %>%
      select(datetime, elapsed, condition, well_id, replicate_id, receptor, treatment, value)
  } else {
    parsed <- guess_receptor_treatment(unique(long$condition))
    
    long %>%
      left_join(parsed, by = c("condition" = "condition_raw")) %>%
      mutate(
        well_id = toupper(extract_well_id(condition)),
        replicate_id = "",
        receptor = if_else(is.na(as.character(receptor)), "none", as.character(receptor)),
        treatment = if_else(is.na(as.character(treatment)), "VEH", as.character(treatment))
      ) %>%
      select(datetime, elapsed, condition, well_id, replicate_id, receptor, treatment, value)
  }
}

# ---------------------------
# Plate map helpers
# ---------------------------
read_plate_map <- function(path) {
  pm <- read_csv(path, show_col_types = FALSE) %>% rename_with(tolower)
  
  required_cols <- c("well", "hormone", "receptor")
  missing <- setdiff(required_cols, names(pm))
  if (length(missing) > 0) {
    stop("Plate map is missing required columns: ", paste(missing, collapse = ", "))
  }
  
  pm %>%
    mutate(
      well = toupper(str_squish(as.character(well))),
      hormone = as.character(hormone),
      receptor = as.character(receptor),
      passage = if ("passage" %in% names(.)) as.character(passage) else NA_character_,
      expt = if ("expt" %in% names(.)) as.character(expt) else NA_character_,
      cell_line = if ("cell_line" %in% names(.)) as.character(cell_line) else NA_character_,
      receptor_pm = purrr::map_chr(receptor, canonicalize_receptor_edit),
      treatment_pm = purrr::map_chr(hormone, canonicalize_treatment_edit),
      passage_pm = if_else(
        is.na(passage) | str_squish(passage) == "",
        NA_character_,
        paste0("Passage_", passage)
      )
    ) %>%
    select(well, cell_line, expt, receptor_pm, treatment_pm, passage_pm)
}

infer_plate_size <- function(wells) {
  wells <- unique(na.omit(wells[wells != ""]))
  if (length(wells) == 0) return(NA_character_)
  
  rows <- str_extract(wells, "^[A-Z]")
  cols <- suppressWarnings(as.integer(str_extract(wells, "[0-9]+$")))
  
  n_rows <- length(unique(rows))
  n_cols <- max(cols, na.rm = TRUE)
  n_total <- n_rows * n_cols
  
  case_when(
    n_rows <= 2 && n_cols <= 3 ~ "6-well",
    n_rows <= 3 && n_cols <= 4 ~ "12-well",
    n_rows <= 4 && n_cols <= 6 ~ "24-well",
    n_rows <= 6 && n_cols <= 8 ~ "48-well",
    n_rows <= 8 && n_cols <= 12 ~ "96-well",
    TRUE ~ paste0(n_total, "-well (inferred)")
  )
}

make_plate_preview_tables <- function(df) {
  df <- df %>%
    mutate(
      well_id = toupper(str_squish(as.character(well_id))),
      plate_id = as.character(passage),
      label = if_else(
        well_id == "" | is.na(well_id),
        NA_character_,
        paste0(
          as.character(receptor), "\n",
          as.character(treatment), "\n",
          as.character(passage)
        )
      )
    ) %>%
    filter(!is.na(well_id), well_id != "") %>%
    distinct(plate_id, well_id, label)
  
  if (nrow(df) == 0) return(list())
  
  split(df, df$plate_id) |>
    purrr::imap(function(d, plate_name) {
      rows <- sort(unique(str_extract(d$well_id, "^[A-Z]")))
      cols <- sort(unique(suppressWarnings(as.integer(str_extract(d$well_id, "[0-9]+$")))))
      
      grid <- expand_grid(row = rows, col = cols) %>%
        mutate(well_id = paste0(row, col)) %>%
        left_join(d %>% select(well_id, label), by = "well_id") %>%
        mutate(label = replace_na(label, "")) %>%
        select(-well_id) %>%
        pivot_wider(names_from = col, values_from = label)
      
      list(
        plate_id = plate_name,
        table = grid
      )
    })
}

# ---------------------------
# Spike masking
# ---------------------------
mask_spikes_cooks <- function(df, value_col = "value_norm", x_col = "elapsed") {
  y <- df[[value_col]]
  x <- df[[x_col]]
  ok <- is.finite(y) & is.finite(x)
  
  df$spike_flag <- FALSE
  df$value_masked <- y
  
  if (sum(ok) < 4) return(df)
  
  fit <- try(lm(y[ok] ~ x[ok]), silent = TRUE)
  if (inherits(fit, "try-error")) return(df)
  
  cooks <- cooks.distance(fit)
  cutoff <- 4 / sum(ok)
  
  spike_idx_local <- which(cooks > cutoff)
  if (length(spike_idx_local) == 0) return(df)
  
  spike_idx_global <- which(ok)[spike_idx_local]
  df$spike_flag[spike_idx_global] <- TRUE
  df$value_masked[spike_idx_global] <- NA_real_
  df
}

# ---------------------------
# Plot/stat helpers
# ---------------------------
preview_tabular_file <- function(path, n = 20) {
  lines <- readLines(path, warn = FALSE)
  header_i <- find_data_header_row(lines)
  if (is.na(header_i)) {
    dat <- read_delim(path, delim = "\t", show_col_types = FALSE, n_max = n)
  } else {
    txt <- paste(lines[header_i:length(lines)], collapse = "\n")
    dat <- read_delim(I(txt), delim = "\t", show_col_types = FALSE, n_max = n)
  }
  dat
}

ols_control_adjust <- function(df, sig_col, ctl_col) {
  x <- df[[ctl_col]]
  y <- df[[sig_col]]
  ok <- is.finite(x) & is.finite(y)
  
  if (sum(ok) < 3) {
    df$value_norm <- NA_real_
    return(df)
  }
  
  fit <- lm(y[ok] ~ x[ok])
  b <- unname(coef(fit)[2])
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

classify_treatment_group <- function(treatment) {
  trt <- toupper(as.character(treatment))
  has_e2   <- str_detect(trt, "\\bE2\\b")
  has_p4   <- str_detect(trt, "\\bP4\\b")
  has_dht  <- str_detect(trt, "\\bDHT\\b")
  has_4oht <- str_detect(trt, "\\b4-OHT\\b")
  has_ru   <- str_detect(trt, "\\bRU-486\\b")
  has_gc   <- str_detect(trt, "\\bDEX\\b|\\bCORT\\b|\\bGLUCO\\b")
  
  ligands <- c(
    if (has_e2) "E2",
    if (has_p4) "P4",
    if (has_dht) "DHT",
    if (has_4oht) "4-OHT",
    if (has_ru) "RU-486",
    if (has_gc) "Glucocorticoid"
  )
  
  if (length(ligands) == 0) return(as.character(treatment))
  paste(ligands, collapse = " + ")
}

treatment_levels_master <- c(
  "VEH","E2","P4","DHT","4-OHT","RU-486","Glucocorticoid",
  "E2 + P4","E2 + DHT","E2 + 4-OHT","E2 + RU-486","P4 + DHT","P4 + 4-OHT","P4 + RU-486",
  "DHT + 4-OHT","DHT + RU-486","4-OHT + RU-486"
)

treatment_color_values <- c(
  "VEH" = "#000000",
  "E2" = "#FB0280",
  "P4" = "#FD8008",
  "DHT" = "#0F80FF",
  "4-OHT" = "#7A3CFF",
  "RU-486" = "#9E9E9E",
  "Glucocorticoid" = "#00A878",
  "E2 + P4" = "#C23B8E",
  "E2 + DHT" = "#8A4DFF",
  "E2 + 4-OHT" = "#B04DFF",
  "E2 + RU-486" = "#B85A8A",
  "P4 + DHT" = "#7F9CFF",
  "P4 + 4-OHT" = "#C06A88",
  "P4 + RU-486" = "#C08A60",
  "DHT + 4-OHT" = "#4F5BFF",
  "DHT + RU-486" = "#6B8BB8",
  "4-OHT + RU-486" = "#8E6BA8"
)

compute_auc_export_dims <- function(df) {
  n_receptors <- dplyr::n_distinct(df$receptor)
  n_treatments <- dplyr::n_distinct(df$treatment_group)
  
  width_mm <- 70 + 12 * n_receptors + 6 * n_treatments
  width_mm <- max(89, min(width_mm, 240))
  height_mm <- 70
  
  list(width_mm = width_mm, height_mm = height_mm)
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  titlePanel("Incucyte Multi-File Import + Channel Normalization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload Incucyte export files (.txt/.tsv/.csv)", multiple = TRUE, accept = c(".txt", ".tsv", ".csv")),
      actionButton("clear_files", "Clear uploaded files", class = "btn-warning"),
      br(), br(),
      checkboxInput("drop_stderr", "Drop '(Std Err ...)' columns", value = TRUE),
      uiOutput("channel_map_ui"),
      hr(),
      uiOutput("norm_ui"),
      actionButton("run", "Import + Process", class = "btn-primary"),
      hr(),
      h4("Downloads"),
      downloadButton("download_prism", "Prism AUC export (csv)"),
      downloadButton("download_timecourse", "Timecourse data (csv)"),
      downloadButton("download_stats_csv", "Export all stats CSV"),
      downloadButton("download_editor", "Editor table (csv)"),
      downloadButton("download_norm", "Normalized by Passage (csv)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          "Plot",
          plotOutput("plot", height = 420),
          br(),
          fluidRow(
            column(8, plotOutput("auc_plot", height = 340)),
            column(
              4,
              br(),
              downloadButton("download_auc_plot_png", "Export AUC plot PNG"),
              br(), br(),
              downloadButton("download_auc_plot_svg", "Export AUC plot SVG")
            )
          ),
          br(),
          fluidRow(
            column(4, uiOutput("plot_time_ui")),
            column(4, uiOutput("plot_receptor_ui")),
            column(4, uiOutput("plot_treatment_ui"))
          )
        ),
        tabPanel(
          "Statistics",
          br(),
          h4("Two-way ANOVA on AUC data"),
          verbatimTextOutput("auc_anova_out"),
          br(),
          h5("Tukey post-hoc: treatments within each receptor"),
          tableOutput("auc_tukey_treat_within_receptor"),
          br(),
          h5("Tukey post-hoc: receptors within each treatment"),
          tableOutput("auc_tukey_receptor_within_treatment"),
          br(),
          hr(),
          h4("3-way repeated-measures ANOVA on time data"),
          verbatimTextOutput("time_anova_out")
        ),
        tabPanel(
          "Factor editor",
          br(),
          fluidRow(
            column(3, actionButton("exclude_selected", "Exclude selected", class = "btn-warning")),
            column(3, actionButton("include_selected", "Include selected")),
            column(3, actionButton("reset_editor", "Reset edits")),
            column(3, actionButton("apply_editor", "Apply factor edits", class = "btn-primary"))
          ),
          br(),
          tags$p("Edit passage, receptor, treatment, or exclude status. Select rows and use exclude/include buttons, or edit the exclude column directly."),
          rHandsontableOutput("editor_table", height = "560px")
        ),
        tabPanel(
          "Plate map / check",
          br(),
          fileInput("platemap", "Upload plate map (.csv)", multiple = FALSE, accept = c(".csv")),
          br(),
          h4("Plate map preview"),
          tableOutput("preview_platemap"),
          br(),
          h4("Plate compatibility summary"),
          verbatimTextOutput("plate_check_summary"),
          br(),
          h4("Plate layout previews"),
          uiOutput("plate_check_layout")
        ),
        tabPanel(
          "File preview",
          br(),
          DTOutput("preview_files_dt")
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  
  uploaded_files_rv <- reactiveVal(NULL)
  
  observeEvent(input$files, {
    req(input$files)
    
    new_files <- tibble(
      file = input$files$name,
      path = input$files$datapath
    )
    
    old_files <- uploaded_files_rv()
    
    if (is.null(old_files)) {
      uploaded_files_rv(new_files)
    } else {
      combined <- bind_rows(old_files, new_files) %>%
        distinct(file, path, .keep_all = TRUE)
      uploaded_files_rv(combined)
    }
  })
  
  observeEvent(input$clear_files, {
    uploaded_files_rv(NULL)
  })
  
  output$channel_map_ui <- renderUI({
    req(uploaded_files_rv())
    fns <- uploaded_files_rv()$file
    
    tagList(
      h4("Assign a fluorescence channel to each file"),
      tags$p(tags$small("You can upload files in batches; files are retained until cleared or the app reloads.")),
      lapply(seq_along(fns), function(i) {
        fname <- fns[i]
        default <- if (str_detect(tolower(fname), "red")) "NIR" else
          if (str_detect(tolower(fname), "nir")) "NIR" else
            if (str_detect(tolower(fname), "gfp|green")) "GFP" else
              if (str_detect(tolower(fname), "orange")) "Orange" else
                "Other"
        
        fluidRow(
          column(8, tags$small(fname)),
          column(4, selectInput(paste0("chan_", i), NULL, c("GFP", "NIR", "Orange", "Red", "Other"), selected = default))
        )
      })
    )
  })
  
  channel_map <- reactive({
    req(uploaded_files_rv())
    files_df <- uploaded_files_rv()
    
    tibble(
      file = files_df$file,
      path = files_df$path,
      channel = purrr::map_chr(seq_len(nrow(files_df)), function(i) {
        fname <- files_df$file[i]
        val <- input[[paste0("chan_", i)]]
        
        if (is.null(val) || is.na(val) || val == "") {
          if (str_detect(tolower(fname), "red")) return("NIR")
          if (str_detect(tolower(fname), "nir")) return("NIR")
          if (str_detect(tolower(fname), "gfp|green")) return("GFP")
          if (str_detect(tolower(fname), "orange")) return("Orange")
          return("Other")
        }
        val
      })
    )
  })
  
  plate_map_tbl <- reactive({
    req(input$platemap)
    read_plate_map(input$platemap$datapath)
  })
  
  output$preview_platemap <- renderTable({
    if (is.null(input$platemap)) return(NULL)
    plate_map_tbl()
  }, striped = TRUE)
  
  output$preview_files_dt <- renderDT({
    req(uploaded_files_rv())
    previews <- purrr::imap(uploaded_files_rv()$path, function(path, i) {
      dat <- preview_tabular_file(path, n = 20)
      dat <- dat %>% mutate(`..file` = uploaded_files_rv()$file[i], .before = 1)
      dat
    })
    bind_rows(previews)
  }, options = list(pageLength = 20, scrollX = TRUE))
  
  output$norm_ui <- renderUI({
    req(channel_map())
    chans <- sort(unique(channel_map()$channel))
    if (length(chans) == 0) return(NULL)
    
    tagList(
      h4("Normalization"),
      selectInput("signal_channel", "Signal channel", choices = chans, selected = chans[1]),
      selectInput("control_channel", "Control channel", choices = chans, selected = if ("NIR" %in% chans) "NIR" else chans[min(2, length(chans))]),
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
      checkboxInput("baseline_norm", "Also baseline-normalize within Passage+Factor combination", value = FALSE),
      checkboxInput("mask_spikes", "Mask spikes (Cook's distance > 4/n)", value = FALSE)
    )
  })
  
  raw_long_auto <- eventReactive(input$run, {
    cm <- channel_map()
    
    purrr::pmap_dfr(cm, function(file, path, channel) {
      meta <- read_incucyte_header_meta(path)
      dat0 <- read_incucyte_long(path, drop_stderr = isTRUE(input$drop_stderr))
      
      default_passage <- if (!is.na(meta$passage[[1]]) && meta$passage[[1]] != "") {
        paste0("Passage_", meta$passage[[1]])
      } else {
        "Passage_NA"
      }
      
      dat0 <- dat0 %>%
        mutate(
          file = file,
          channel = channel,
          well_id = toupper(str_squish(well_id)),
          passage = if_else(replicate_id != "", paste0(default_passage, "_rep", replicate_id), default_passage),
          vessel_name = meta$vessel_name[[1]] %||% NA_character_,
          metric      = meta$metric[[1]] %||% NA_character_,
          cell_type   = meta$cell_type[[1]] %||% NA_character_,
          analysis    = meta$analysis[[1]] %||% NA_character_,
          condition_id = paste(file, condition, sep = " || "),
          receptor = purrr::map_chr(receptor, canonicalize_receptor_combo),
          treatment = purrr::map_chr(treatment, canonicalize_treatment_combo)
        )
      
      if (!is.null(input$platemap)) {
        dat0 <- dat0 %>%
          left_join(plate_map_tbl(), by = c("well_id" = "well")) %>%
          mutate(
            receptor = coalesce(receptor_pm, receptor),
            treatment = coalesce(treatment_pm, treatment),
            passage = coalesce(passage_pm, passage)
          ) %>%
          select(-cell_line, -expt, -receptor_pm, -treatment_pm, -passage_pm)
      }
      
      dat0 %>%
        mutate(factor_key = make_factor_key(receptor, treatment)) %>%
        select(
          condition_id, file, well_id, replicate_id, channel, passage,
          vessel_name, metric, cell_type, analysis,
          datetime, elapsed, condition, receptor, treatment, factor_key, value
        )
    })
  }, ignoreInit = TRUE)
  
  output$plate_check_summary <- renderPrint({
    req(raw_long_auto())
    wells_in_data <- raw_long_auto() %>% pull(well_id) %>% unique()
    plate_size_data <- infer_plate_size(wells_in_data)
    
    if (is.null(input$platemap)) {
      cat("Detected data plate size:", plate_size_data, "\n")
      cat("No plate map uploaded.\n")
    } else {
      pm <- plate_map_tbl()
      wells_pm <- pm$well
      plate_size_pm <- infer_plate_size(wells_pm)
      
      matched <- intersect(wells_in_data[wells_in_data != ""], wells_pm)
      unmatched_data <- setdiff(wells_in_data[wells_in_data != ""], wells_pm)
      unmatched_pm <- setdiff(wells_pm, wells_in_data[wells_in_data != ""])
      
      cat("Detected data plate size:", plate_size_data, "\n")
      cat("Detected plate map size:", plate_size_pm, "\n")
      cat("Matched wells:", length(matched), "\n")
      cat("Unmatched data wells:", length(unmatched_data), "\n")
      cat("Unmatched plate map wells:", length(unmatched_pm), "\n")
      if (length(unmatched_data) > 0) cat("Data-only wells:", paste(sort(unmatched_data), collapse = ", "), "\n")
      if (length(unmatched_pm) > 0) cat("Plate-map-only wells:", paste(sort(unmatched_pm), collapse = ", "), "\n")
    }
  })
  
  editor_default <- reactive({
    req(raw_long_auto())
    
    rep_counts <- raw_long_auto() %>%
      distinct(condition, file) %>%
      count(condition, name = "n_reps")
    
    raw_long_auto() %>%
      distinct(condition_id, file, well_id, replicate_id, condition, passage, receptor, treatment, factor_key) %>%
      left_join(rep_counts, by = "condition") %>%
      mutate(
        exclude = FALSE,
        well_id = as.character(well_id),
        replicate_id = as.character(replicate_id),
        passage = as.character(passage),
        receptor = as.character(receptor),
        treatment = as.character(treatment),
        factor_key = as.character(factor_key),
        n_reps = as.integer(n_reps)
      ) %>%
      arrange(file, well_id, replicate_id, condition)
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
      select(exclude, file, well_id, replicate_id, condition, n_reps, passage, receptor, treatment)
    
    rhandsontable(df, rowHeaders = NULL, stretchH = "all", height = 540, selectCallback = TRUE) %>%
      hot_col("exclude", type = "checkbox") %>%
      hot_col("file", readOnly = TRUE) %>%
      hot_col("well_id", readOnly = TRUE) %>%
      hot_col("replicate_id", readOnly = TRUE) %>%
      hot_col("condition", readOnly = TRUE) %>%
      hot_col("n_reps", readOnly = TRUE) %>%
      hot_table(highlightCol = TRUE, highlightRow = TRUE, columnSorting = TRUE, manualColumnMove = TRUE)
  })
  
  observe({
    if (!is.null(input$editor_table)) {
      tbl <- hot_to_r(input$editor_table)
      if (!is.null(tbl) && !is.null(editor_rv())) {
        tbl <- as_tibble(tbl) %>%
          mutate(
            exclude = as.logical(exclude),
            file = as.character(file),
            well_id = as.character(well_id),
            replicate_id = as.character(replicate_id),
            condition = as.character(condition),
            n_reps = suppressWarnings(as.integer(n_reps)),
            passage = as.character(passage),
            receptor = as.character(receptor),
            treatment = as.character(treatment)
          )
        
        key_map <- editor_rv() %>%
          transmute(condition_id, file, well_id, replicate_id, condition) %>%
          distinct()
        
        updated <- tbl %>%
          left_join(key_map, by = c("file", "well_id", "replicate_id", "condition")) %>%
          select(condition_id, everything())
        
        editor_rv(updated)
      }
    }
  })
  
  observeEvent(input$exclude_selected, {
    req(editor_rv(), input$editor_table_select$select$r)
    rows <- input$editor_table_select$select$r
    df <- editor_rv()
    rows <- rows[rows >= 1 & rows <= nrow(df)]
    if (length(rows) > 0) {
      df$exclude[rows] <- TRUE
      editor_rv(df)
    }
  })
  
  observeEvent(input$include_selected, {
    req(editor_rv(), input$editor_table_select$select$r)
    rows <- input$editor_table_select$select$r
    df <- editor_rv()
    rows <- rows[rows >= 1 & rows <= nrow(df)]
    if (length(rows) > 0) {
      df$exclude[rows] <- FALSE
      editor_rv(df)
    }
  })
  
  edited_map <- reactiveVal(NULL)
  
  observeEvent(input$apply_editor, {
    req(editor_rv())
    
    df_edit <- editor_rv() %>%
      mutate(
        exclude = if_else(is.na(exclude), FALSE, as.logical(exclude)),
        well_id = as.character(well_id),
        replicate_id = as.character(replicate_id),
        passage = purrr::map_chr(passage, canonicalize_passage_edit),
        receptor = purrr::map_chr(receptor, canonicalize_receptor_edit),
        treatment = purrr::map_chr(treatment, canonicalize_treatment_edit),
        factor_key = make_factor_key(receptor, treatment),
        receptor = factor(receptor),
        treatment = factor(treatment),
        passage = factor(passage)
      ) %>%
      arrange(file, well_id, replicate_id, condition)
    
    edited_map(df_edit)
  }, ignoreInit = TRUE)
  
  current_editor_map <- reactive({
    if (!is.null(edited_map())) {
      edited_map()
    } else {
      req(editor_default())
      editor_default() %>%
        mutate(
          exclude = if_else(is.na(exclude), FALSE, as.logical(exclude)),
          well_id = as.character(well_id),
          replicate_id = as.character(replicate_id),
          passage = factor(purrr::map_chr(passage, canonicalize_passage_edit)),
          receptor = factor(purrr::map_chr(receptor, canonicalize_receptor_edit)),
          treatment = factor(purrr::map_chr(treatment, canonicalize_treatment_edit)),
          factor_key = make_factor_key(receptor, treatment)
        )
    }
  })
  
  output$plate_check_layout <- renderUI({
    req(current_editor_map())
    
    plate_tables <- make_plate_preview_tables(current_editor_map())
    
    if (length(plate_tables) == 0) {
      return(tags$p("No well-based plate layout available."))
    }
    
    tagList(
      purrr::map(plate_tables, function(x) {
        tagList(
          tags$h4(paste("Plate:", x$plate_id)),
          tableOutput(outputId = paste0("plate_layout_", make.names(x$plate_id))),
          tags$br()
        )
      })
    )
  })
  
  observe({
    req(current_editor_map())
    plate_tables <- make_plate_preview_tables(current_editor_map())
    
    purrr::walk(plate_tables, function(x) {
      local({
        plate_name <- x$plate_id
        plate_tbl <- x$table
        output_id <- paste0("plate_layout_", make.names(plate_name))
        
        output[[output_id]] <- renderTable({
          plate_tbl
        }, striped = TRUE, bordered = TRUE, spacing = "xs")
      })
    })
  })
  
  raw_long <- reactive({
    req(raw_long_auto(), current_editor_map())
    
    raw_long_auto() %>%
      select(-passage, -receptor, -treatment, -factor_key) %>%
      left_join(
        current_editor_map() %>% select(condition_id, exclude, passage, receptor, treatment, factor_key),
        by = "condition_id"
      ) %>%
      mutate(
        exclude = if_else(is.na(exclude), FALSE, as.logical(exclude)),
        passage = factor(as.character(passage)),
        receptor = factor(purrr::map_chr(as.character(receptor), canonicalize_receptor_combo)),
        treatment = factor(as.character(treatment)),
        factor_key = make_factor_key(receptor, treatment)
      ) %>%
      filter(!exclude) %>%
      relocate(exclude, well_id, replicate_id, passage, receptor, treatment, factor_key, .after = condition)
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
      out <- w %>% mutate(value_norm = NA_real_)
    } else {
      out <- w
      
      if (method %in% c("none", "ratio", "log2ratio")) {
        out <- out %>%
          mutate(
            value_norm = case_when(
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
    }
    
    if (isTRUE(input$mask_spikes)) {
      out <- out %>%
        group_by(passage, factor_key) %>%
        group_modify(~ mask_spikes_cooks(.x, value_col = "value_norm", x_col = "elapsed")) %>%
        ungroup() %>%
        mutate(value_norm = value_masked) %>%
        select(-value_masked)
    } else {
      out <- out %>% mutate(spike_flag = FALSE)
    }
    
    out
  })
  
  stats_long <- reactive({
    req(normalized_passage())
    normalized_passage() %>%
      select(passage, elapsed, receptor, treatment, factor_key, value_norm, spike_flag) %>%
      arrange(receptor, treatment, passage, elapsed)
  })
  
  output$plot_time_ui <- renderUI({
    req(stats_long())
    df <- stats_long() %>% filter(is.finite(elapsed))
    if (nrow(df) == 0) return(NULL)
    
    rng <- range(df$elapsed, na.rm = TRUE)
    minv <- floor(rng[1])
    maxv <- ceiling(rng[2])
    
    sliderInput(
      "plot_time_range",
      "Elapsed time range",
      min = minv,
      max = maxv,
      value = c(minv, maxv),
      step = 1
    )
  })
  
  output$plot_receptor_ui <- renderUI({
    req(stats_long())
    levs <- stats_long() %>% pull(receptor) %>% as.character() %>% unique() %>% sort()
    checkboxGroupInput("plot_receptors", "Receptors to show", choices = levs, selected = levs)
  })
  
  output$plot_treatment_ui <- renderUI({
    req(stats_long())
    levs <- stats_long() %>% pull(treatment) %>% as.character() %>% unique() %>% sort()
    checkboxGroupInput("plot_treatments", "Treatments to show", choices = levs, selected = levs)
  })
  
  filtered_stats_long <- reactive({
    req(stats_long())
    df <- stats_long()
    
    if (!is.null(input$plot_time_range)) {
      df <- df %>% filter(elapsed >= input$plot_time_range[1], elapsed <= input$plot_time_range[2])
    }
    if (!is.null(input$plot_receptors) && length(input$plot_receptors) > 0) {
      df <- df %>% filter(as.character(receptor) %in% input$plot_receptors)
    }
    if (!is.null(input$plot_treatments) && length(input$plot_treatments) > 0) {
      df <- df %>% filter(as.character(treatment) %in% input$plot_treatments)
    }
    
    df
  })
  
  auc_preview <- reactive({
    req(filtered_stats_long())
    filtered_stats_long() %>%
      group_by(receptor, treatment, passage) %>%
      summarize(auc = auc_trapz(elapsed, value_norm), .groups = "drop")
  })
  
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
  
  auc_plot_obj <- reactive({
    req(auc_preview())
    df <- auc_preview()
    if (nrow(df) == 0 || all(!is.finite(df$auc))) return(NULL)
    
    df <- df %>%
      mutate(
        treatment_group = purrr::map_chr(treatment, classify_treatment_group),
        treatment_group = factor(
          treatment_group,
          levels = unique(c(
            treatment_levels_master,
            sort(setdiff(unique(treatment_group), treatment_levels_master))
          ))
        )
      )
    
    missing_levels <- setdiff(levels(df$treatment_group), names(treatment_color_values))
    color_values <- treatment_color_values
    if (length(missing_levels) > 0) {
      extra_cols <- rep("#666666", length(missing_levels))
      names(extra_cols) <- missing_levels
      color_values <- c(color_values, extra_cols)
    }
    
    dodge <- position_dodge(width = 0.8)
    
    ggplot(df, aes(x = receptor, y = auc, color = treatment_group)) +
      geom_point(position = dodge, size = 2.8, alpha = 0.5, stroke = 0) +
      scale_color_manual(values = color_values, drop = TRUE) +
      labs(
        x = "Receptor",
        y = "AUC",
        color = "Treatment",
        title = "Combined AUC preview for current filter window"
      ) +
      theme_classic(base_size = 8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.4),
        axis.ticks = element_line(linewidth = 0.4)
      )
  })
  
  output$auc_plot <- renderPlot({
    p <- auc_plot_obj()
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "No AUC values available for the current filters.")
      return()
    }
    p
  })
  
  # ---------------------------
  # Statistics
  # ---------------------------
  time_stats_df <- reactive({
    req(filtered_stats_long())
    
    filtered_stats_long() %>%
      mutate(
        receptor = factor(receptor),
        treatment = factor(treatment),
        elapsed_f = factor(elapsed),
        subject_id = factor(passage)
      ) %>%
      filter(is.finite(value_norm))
  })
  
  time_rm_fit <- reactive({
    df <- time_stats_df()
    
    validate(
      need(nrow(df) > 0, "No time-course data available for repeated-measures ANOVA."),
      need(n_distinct(df$receptor) >= 2, "Need at least 2 receptor levels for time ANOVA."),
      need(n_distinct(df$treatment) >= 2, "Need at least 2 treatment levels for time ANOVA."),
      need(n_distinct(df$elapsed_f) >= 2, "Need at least 2 elapsed timepoints for time ANOVA."),
      need(n_distinct(df$subject_id) >= 2, "Need at least 2 repeated subjects/replicates for time ANOVA.")
    )
    
    aov(value_norm ~ receptor * treatment * elapsed_f + Error(subject_id / elapsed_f), data = df)
  })
  
  auc_stats_df <- reactive({
    req(auc_preview())
    auc_preview() %>%
      mutate(
        receptor = factor(receptor),
        treatment = factor(treatment)
      ) %>%
      filter(is.finite(auc))
  })
  
  auc_anova_fit <- reactive({
    df <- auc_stats_df()
    
    validate(
      need(nrow(df) > 0, "No AUC data available for ANOVA."),
      need(n_distinct(df$receptor) >= 2, "Need at least 2 receptor levels for AUC ANOVA."),
      need(n_distinct(df$treatment) >= 2, "Need at least 2 treatment levels for AUC ANOVA.")
    )
    
    aov(auc ~ receptor * treatment, data = df)
  })
  
  auc_tukey_treatment_within_receptor_df <- reactive({
    df <- auc_stats_df()
    
    out <- split(df, df$receptor) |>
      purrr::imap_dfr(function(d, rec_lab) {
        if (n_distinct(d$treatment) < 2) return(NULL)
        fit <- try(aov(auc ~ treatment, data = d), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        tk <- try(TukeyHSD(fit, "treatment"), silent = TRUE)
        if (inherits(tk, "try-error")) return(NULL)
        
        as.data.frame(tk$treatment) %>%
          tibble::rownames_to_column("comparison") %>%
          mutate(receptor = as.character(rec_lab), .before = 1)
      })
    
    if (nrow(out) == 0) {
      tibble(receptor = character(), comparison = character(), diff = numeric(), lwr = numeric(), upr = numeric(), `p adj` = numeric())
    } else out
  })
  
  auc_tukey_receptor_within_treatment_df <- reactive({
    df <- auc_stats_df()
    
    out <- split(df, df$treatment) |>
      purrr::imap_dfr(function(d, trt_lab) {
        if (n_distinct(d$receptor) < 2) return(NULL)
        fit <- try(aov(auc ~ receptor, data = d), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        tk <- try(TukeyHSD(fit, "receptor"), silent = TRUE)
        if (inherits(tk, "try-error")) return(NULL)
        
        as.data.frame(tk$receptor) %>%
          tibble::rownames_to_column("comparison") %>%
          mutate(treatment = as.character(trt_lab), .before = 1)
      })
    
    if (nrow(out) == 0) {
      tibble(treatment = character(), comparison = character(), diff = numeric(), lwr = numeric(), upr = numeric(), `p adj` = numeric())
    } else out
  })
  
  output$auc_anova_out <- renderPrint({
    summary(auc_anova_fit())
  })
  
  output$auc_tukey_treat_within_receptor <- renderTable({
    auc_tukey_treatment_within_receptor_df()
  }, striped = TRUE)
  
  output$auc_tukey_receptor_within_treatment <- renderTable({
    auc_tukey_receptor_within_treatment_df()
  }, striped = TRUE)
  
  output$time_anova_out <- renderPrint({
    summary(time_rm_fit())
  })
  
  # ---------------------------
  # Downloads
  # ---------------------------
  output$download_editor <- downloadHandler(
    filename = function() paste0("incucyte_editor_table_", Sys.Date(), ".csv"),
    content = function(file) write_csv(current_editor_map(), file)
  )
  
  output$download_norm <- downloadHandler(
    filename = function() paste0("incucyte_normalized_by_passage_", Sys.Date(), ".csv"),
    content = function(file) write_csv(normalized_passage(), file)
  )
  
  output$download_prism <- downloadHandler(
    filename = function() paste0("incucyte_prism_auc_", Sys.Date(), ".csv"),
    content = function(file) {
      auc_df <- auc_preview() %>% arrange(receptor, treatment, passage)
      
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
        pivot_wider(names_from = col_key, values_from = auc) %>%
        select(receptor, any_of(col_template))
      
      value_cols <- names(wide)[-1]
      treatment_header <- c("", str_replace(value_cols, "__rep\\d+$", ""))
      rep_header <- c("receptor", str_extract(value_cols, "rep\\d+$"))
      
      out <- rbind(treatment_header, rep_header, as.matrix(wide))
      
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
  
  output$download_timecourse <- downloadHandler(
    filename = function() paste0("incucyte_timecourse_", Sys.Date(), ".csv"),
    content = function(file) {
      
      req(filtered_stats_long())
      
      df <- filtered_stats_long() %>%
        mutate(
          cell_line = as.factor(cell_type %||% "unknown"),
          hormone   = factor(as.character(treatment), ordered = TRUE),
          receptor  = factor(as.character(receptor)),
          expt      = factor(str_extract(as.character(file), "[0-9]+") %||% "1"),
          passage   = factor(str_replace(as.character(passage), "Passage_", "p")),
          elapsed_hour = elapsed,
          id = factor(paste0(expt, "_", passage))
        ) %>%
        group_by(cell_line, hormone, receptor, expt, passage, elapsed_hour, id) %>%
        summarise(
          mean_count   = mean(value_norm, na.rm = TRUE),
          median_count = median(value_norm, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        arrange(cell_line, hormone, receptor, expt, passage, elapsed_hour)
      
      write_csv(df, file)
    }
  )
  
  output$download_auc_plot_png <- downloadHandler(
    filename = function() paste0("incucyte_auc_plot_", Sys.Date(), ".png"),
    content = function(file) {
      p <- auc_plot_obj()
      if (is.null(p)) {
        png(file, width = 1800, height = 1200, res = 300)
        plot.new()
        text(0.5, 0.5, "No AUC values available for the current filters.")
        dev.off()
        return()
      }
      
      df <- auc_preview() %>%
        mutate(treatment_group = purrr::map_chr(treatment, classify_treatment_group))
      dims <- compute_auc_export_dims(df)
      
      ggsave(
        filename = file,
        plot = p,
        device = "png",
        width = dims$width_mm,
        height = dims$height_mm,
        units = "mm",
        dpi = 600,
        bg = "white"
      )
    }
  )
  
  output$download_auc_plot_svg <- downloadHandler(
    filename = function() paste0("incucyte_auc_plot_", Sys.Date(), ".svg"),
    content = function(file) {
      p <- auc_plot_obj()
      if (is.null(p)) {
        svg(filename = file, width = 3.5, height = 2.75, bg = "white")
        plot.new()
        text(0.5, 0.5, "No AUC values available for the current filters.")
        dev.off()
        return()
      }
      
      df <- auc_preview() %>%
        mutate(treatment_group = purrr::map_chr(treatment, classify_treatment_group))
      dims <- compute_auc_export_dims(df)
      
      ggsave(
        filename = file,
        plot = p,
        device = "svg",
        width = dims$width_mm,
        height = dims$height_mm,
        units = "mm",
        bg = "white"
      )
    }
  )
  
  output$download_stats_csv <- downloadHandler(
    filename = function() paste0("incucyte_stats_", Sys.Date(), ".csv"),
    content = function(file) {
      auc_anova_tbl <- broom::tidy(auc_anova_fit()) %>%
        mutate(table = "auc_anova", .before = 1)
      
      time_anova_tbl <- tryCatch({
        broom::tidy(time_rm_fit()) %>%
          mutate(table = "time_anova", .before = 1)
      }, error = function(e) {
        tibble(table = "time_anova", term = NA_character_, statistic = NA_real_, p.value = NA_real_)
      })
      
      auc_treat_tbl <- auc_tukey_treatment_within_receptor_df() %>%
        mutate(table = "auc_tukey_treatment_within_receptor", .before = 1)
      
      auc_receptor_tbl <- auc_tukey_receptor_within_treatment_df() %>%
        mutate(table = "auc_tukey_receptor_within_treatment", .before = 1)
      
      bind_rows(
        suppressWarnings(auc_anova_tbl),
        suppressWarnings(time_anova_tbl),
        suppressWarnings(auc_treat_tbl),
        suppressWarnings(auc_receptor_tbl)
      ) %>%
        write_csv(file)
    }
  )
}

shinyApp(ui, server)