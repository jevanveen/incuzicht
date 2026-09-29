# app.R
# Incucyte Multi-File Import + Channel Normalization
# HITL factor-map editor, buffered edits, plate map preview, plotting, AUC, exports

library(shiny)
library(tidyverse)
library(readr)
library(stringr)
library(DT)

`%||%` <- function(x, y) if (is.null(x)) y else x

# ---------------------------
# Low-level helpers
# ---------------------------
find_data_header_row <- function(lines) {
  idx <- which(str_detect(lines, "^Date\\s*Time\tElapsed\\b"))
  if (length(idx) == 0) return(NA_integer_)
  idx[1]
}

extract_key_value <- function(lines, key) {
  pat <- paste0("^", key, "\\s*:\\s*(.*)$")
  hit <- lines[str_detect(lines, pat)]
  if (length(hit) == 0) return(NA_character_)
  val <- str_match(hit[1], pat)[, 2]
  if (is.na(val)) NA_character_ else str_trim(val)
}

extract_well_id <- function(x) {
  z <- toupper(str_squish(as.character(x)))
  m_paren <- str_match(z, "\\(\\s*([A-H]\\s*0?(?:[1-9]|1[0-2]))\\s*\\)")
  m_bare <- str_match(z, "^([A-H]\\s*0?(?:[1-9]|1[0-2]))$")
  out <- coalesce(m_paren[, 2], m_bare[, 2], "")
  out <- str_replace_all(out, "\\s+", "")
  out <- str_replace(out, "^([A-H])0", "\\1")
  out[is.na(out)] <- ""
  out
}

standardize_passage_label <- function(x) {
  z <- str_trim(as.character(x))
  if (length(z) != 1 || is.na(z) || z == "") return("Passage_NA")
  if (str_detect(z, regex("^Passage_", ignore_case = TRUE))) return(z)
  if (str_detect(z, regex("^p[0-9]+$", ignore_case = TRUE))) return(paste0("Passage_", str_extract(z, "[0-9]+")))
  if (str_detect(z, "^[0-9]+$")) return(paste0("Passage_", z))
  z
}

clean_user_factor <- function(x) {
  str_squish(as.character(x))
}

read_incucyte_file <- function(path, drop_stderr = TRUE) {
  lines <- readLines(path, warn = FALSE)
  header_i <- find_data_header_row(lines)
  
  if (is.na(header_i)) {
    stop("Couldn't find data header row starting with 'Date Time<TAB>Elapsed'.")
  }
  
  meta_lines <- lines[1:(header_i - 1)]
  
  meta <- tibble(
    vessel_name = extract_key_value(meta_lines, "Vessel Name"),
    metric      = extract_key_value(meta_lines, "Metric"),
    cell_type   = extract_key_value(meta_lines, "Cell Type"),
    passage     = extract_key_value(meta_lines, "Passage"),
    notes       = extract_key_value(meta_lines, "Notes"),
    analysis    = extract_key_value(meta_lines, "Analysis")
  )
  
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
    dat <- dat %>% select(-matches("Std Err"))
  }
  
  list(meta = meta, dat = dat)
}

guess_channel <- function(file) {
  case_when(
    str_detect(tolower(file), "red|nir") ~ "NIR",
    str_detect(tolower(file), "gfp|green") ~ "GFP",
    str_detect(tolower(file), "orange") ~ "Orange",
    TRUE ~ "Other"
  )
}

plate_source_key <- function(file, parsed) {
  stem <- tools::file_path_sans_ext(basename(file)) %>%
    str_to_lower() %>%
    str_remove("[ _-]*(green|gfp|red|nir|orange)$")
  wells <- extract_well_id(names(parsed$dat)[-c(1, 2)])
  layout <- if (all(wells != "")) sort(wells) else sort(names(parsed$dat)[-c(1, 2)])
  paste(stem, parsed$meta$vessel_name[[1]], parsed$dat$datetime[1],
        paste(layout, collapse = "|"), sep = " || ")
}

join_plate_map <- function(long, plate_map, plate_ids) {
  if ("plate_id" %in% names(plate_map)) {
    if (any(is.na(plate_map$plate_id) | plate_map$plate_id == "")) {
      stop("Every map row must have a plate_id when that column is supplied.")
    }
    if (any(!plate_map$plate_id %in% plate_ids)) {
      stop("Plate map contains unknown plate identifiers. Check plate assignments or export a matching map.")
    }
    join_keys <- c("plate_id", "well_id" = "well")
    map_keys <- c("plate_id", "well")
  } else {
    if (length(plate_ids) > 1) {
      stop("A map without plate_id can only be used with one plate. Export a multi-plate .zicht map first.")
    }
    join_keys <- c("well_id" = "well")
    map_keys <- "well"
  }
  if (anyDuplicated(plate_map[map_keys])) stop("Plate map contains duplicate plate/well rows.")
  long %>%
    left_join(plate_map, by = join_keys) %>%
    mutate(
      receptor_pm = coalesce(receptor_pm.y, receptor_pm.x),
      treatment_pm = coalesce(treatment_pm.y, treatment_pm.x),
      passage_pm = coalesce(passage_pm.y, passage_pm.x),
      cell_line_pm = coalesce(cell_line, cell_line_pm),
      expt_pm = coalesce(expt, expt_pm),
      passage = coalesce(passage_pm, passage)
    ) %>%
    select(-ends_with("_pm.x"), -ends_with("_pm.y"), -cell_line, -expt)
}

# ---------------------------
# Canonicalization for initial guesses only
# ---------------------------
canonicalize_receptor_one <- function(x) {
  if (length(x) != 1) stop("canonicalize_receptor_one expects scalar input")
  if (is.na(x) || x == "" || tolower(x) == "none") return("none")
  
  parts <- str_split(as.character(x), "\\s*\\+\\s*", simplify = FALSE)[[1]]
  parts <- str_squish(parts)
  parts <- parts[parts != ""]
  parts <- unique(parts)
  
  canon_part <- function(p) {
    z <- p %>%
      str_to_lower() %>%
      str_replace_all("\u03b1", "a") %>%
      str_replace_all("\u03b2", "b") %>%
      str_replace_all("[^a-z0-9_ ]", " ") %>%
      str_squish()
    
    case_when(
      z %in% c("era", "er", "er a", "esr1", "er1", "er_a") ~ "ER_a",
      z %in% c("erb", "er b", "esr2", "er2", "er_b")       ~ "ER_b",
      z %in% c("pgr", "pr")                                 ~ "PR",
      z %in% c("pra", "pr a", "pr_a")                       ~ "PR_a",
      z %in% c("prb", "pr b", "pr_b")                       ~ "PR_b",
      z %in% c("ar")                                        ~ "AR",
      z %in% c("gr")                                        ~ "GR",
      z %in% c("mr")                                        ~ "MR",
      z %in% c("", "none", "na")                            ~ "none",
      TRUE                                                  ~ str_trim(p)
    )
  }
  
  parts <- purrr::map_chr(parts, canon_part)
  parts <- unique(parts[parts != "none" & parts != ""])
  
  if (any(c("PR_a", "PR_b") %in% parts)) {
    parts <- setdiff(parts, "PR")
  }
  
  ord <- c("ER_a", "ER_b", "PR", "PR_a", "PR_b", "AR", "GR", "MR")
  parts <- parts[order(match(parts, ord, nomatch = 999), parts)]
  
  if (length(parts) == 0) "none" else paste(parts, collapse = " + ")
}

canonicalize_receptor_combo <- Vectorize(canonicalize_receptor_one, USE.NAMES = FALSE)

canonicalize_treatment_one <- function(x) {
  z <- str_trim(as.character(x))
  if (length(z) != 1 || is.na(z) || z == "" || toupper(z) == "VEH") return("VEH")
  
  parts <- str_split(z, "\\s*\\+\\s*", simplify = FALSE)[[1]]
  parts <- toupper(str_squish(parts))
  parts <- str_replace_all(parts, "\\bRU\\s*486\\b", "RU-486")
  parts <- str_replace_all(parts, "\\bFULVESTRANT\\b", "FUL")
  parts <- parts[parts != ""]
  parts <- unique(parts)
  
  get_ligand <- function(s) {
    m <- str_match(s, "\\b(E2|P4|DHT|4-OHT|DEX|CORT|RU-486|NET|FUL)\\b")
    lig <- m[, 2]
    ifelse(is.na(lig), "ZZZ", lig)
  }
  
  get_dose <- function(s) {
    m <- str_match(s, "^([0-9]+\\.?[0-9]*)")
    dose <- suppressWarnings(as.numeric(m[, 2]))
    ifelse(is.na(dose), Inf, dose)
  }
  
  ligand_order <- c("E2", "P4", "DHT", "4-OHT", "NET", "RU-486", "FUL", "DEX", "CORT", "ZZZ")
  ord_lig <- match(get_ligand(parts), ligand_order)
  ord_dose <- get_dose(parts)
  
  parts <- parts[order(ord_lig, ord_dose, parts)]
  paste(parts, collapse = " + ")
}

canonicalize_treatment_combo <- Vectorize(canonicalize_treatment_one, USE.NAMES = FALSE)

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

make_guess_key <- function(receptor_guess, treatment_guess, passage_guess) {
  paste(
    normalize_factor_key(receptor_guess),
    normalize_factor_key(treatment_guess),
    normalize_factor_key(passage_guess),
    sep = " || "
  )
}

# ---------------------------
# Parsing
# ---------------------------
clean_condition_header <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("\\([a-h][0-9]{1,2}\\)", " ") %>%
    str_replace_all(",", " + ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*ul\\b", " ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*k\\s*/\\s*well\\b", " ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*mg/ml\\b", " ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*ug/ml\\b", " ") %>%
    str_replace_all("\\b[0-9]+\\.?[0-9]*\\s*µg/ml\\b", " ") %>%
    str_replace_all("([0-9]+\\.?[0-9]*)(pm|nm|um|\u00b5m|mm)\\b", "\\1 \\2") %>%
    str_replace_all("\\bru\\s*486\\b", "ru-486") %>%
    str_replace_all("\\bfulvestrant\\b", "ful") %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s*\\+\\s*", " + ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
}

extract_receptors <- function(x) {
  s <- clean_condition_header(x)
  recs <- c()
  
  if (str_detect(s, "\\berb\\b|\\ber b\\b|\\besr2\\b|\\ber2\\b")) recs <- c(recs, "ER_b")
  
  if (str_detect(s, "\\bera\\b|\\ber a\\b|\\besr1\\b|\\ber1\\b")) {
    recs <- c(recs, "ER_a")
  } else if (str_detect(s, "\\ber\\b") && !str_detect(s, "\\berb\\b|\\ber b\\b")) {
    recs <- c(recs, "ER_a")
  }
  
  if (str_detect(s, "\\bpra\\b|\\bpr a\\b")) recs <- c(recs, "PR_a")
  if (str_detect(s, "\\bprb\\b|\\bpr b\\b")) recs <- c(recs, "PR_b")
  
  if (
    str_detect(s, "\\bpgr\\b") ||
    (str_detect(s, "\\bpr\\b") && !str_detect(s, "\\bpra\\b|\\bpr a\\b|\\bprb\\b|\\bpr b\\b"))
  ) {
    recs <- c(recs, "PR")
  }
  
  if (str_detect(s, "\\bar\\b|\\bandrogen receptor\\b|\\bnr3c4\\b")) recs <- c(recs, "AR")
  if (str_detect(s, "\\bgr\\b|\\bglucocorticoid receptor\\b|\\bnr3c1\\b")) recs <- c(recs, "GR")
  if (str_detect(s, "\\bmr\\b|\\bmineralocorticoid receptor\\b|\\bnr3c2\\b")) recs <- c(recs, "MR")
  
  if (length(recs) == 0) return("none")
  canonicalize_receptor_one(paste(recs, collapse = " + "))
}

extract_treatment <- function(x) {
  s <- clean_condition_header(x)
  
  veh_pats <- c("\\bveh\\b", "\\bvehicle\\b", "\\be2oh\\b", "\\bethanol\\b", "\\betoh\\b", "\\bdmso\\b")
  veh_only <- any(purrr::map_lgl(veh_pats, ~ str_detect(s, .x)))
  
  tokens <- str_split(s, "\\s+", simplify = TRUE)
  tokens <- tokens[tokens != ""]
  
  is_num <- function(z) str_detect(z, "^[0-9]+\\.?[0-9]*$")
  is_unit <- function(z) str_detect(z, regex("^(pm|nm|um|\u00b5m|mm)$", ignore_case = TRUE))
  is_ligand <- function(z) str_detect(z, regex("^(e2|p4|dht|4-oht|dex|cort|ru-486|net|ful)$", ignore_case = TRUE))
  
  norm_unit <- function(z) {
    z <- toupper(z)
    if (z == "\u00b5M") "uM" else z
  }
  
  norm_lig <- function(z) toupper(z)
  
  hits <- character(0)
  used <- rep(FALSE, length(tokens))
  
  i <- 1
  while (i <= length(tokens)) {
    if (used[i]) {
      i <- i + 1
      next
    }
    
    if (
      i + 2 <= length(tokens) &&
      !any(used[i:(i + 2)]) &&
      is_ligand(tokens[i]) &&
      is_num(tokens[i + 1]) &&
      is_unit(tokens[i + 2])
    ) {
      hits <- c(hits, paste0(tokens[i + 1], norm_unit(tokens[i + 2]), " ", norm_lig(tokens[i])))
      used[i:(i + 2)] <- TRUE
      i <- i + 3
      next
    }
    
    if (
      i + 2 <= length(tokens) &&
      !any(used[i:(i + 2)]) &&
      is_num(tokens[i]) &&
      is_unit(tokens[i + 1]) &&
      is_ligand(tokens[i + 2])
    ) {
      hits <- c(hits, paste0(tokens[i], norm_unit(tokens[i + 1]), " ", norm_lig(tokens[i + 2])))
      used[i:(i + 2)] <- TRUE
      i <- i + 3
      next
    }
    
    i <- i + 1
  }
  
  hits <- unique(hits)
  if (length(hits) > 0) return(canonicalize_treatment_one(paste(hits, collapse = " + ")))
  if (veh_only) return("VEH")
  
  s2 <- s %>%
    str_replace_all(
      "\\b(era|er a|erb|er b|er|esr1|esr2|pra|pr a|prb|pr b|pgr|pr|ar|gr|mr|androgen receptor|glucocorticoid receptor|mineralocorticoid receptor|hek293t|293t|noer|reporter)\\b",
      " "
    ) %>%
    str_replace_all("\\s*\\+\\s*", " + ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim()
  
  if (s2 == "" || s2 == "+") return("VEH")
  toupper(s2)
}

extract_passage_from_condition <- function(x) {
  m <- str_match(x, "\\bp([0-9]+)\\b")
  if (!is.na(m[1, 2])) return(paste0("Passage_", m[1, 2]))
  NA_character_
}

map_distinct_chr <- function(values, transform) {
  distinct_values <- unique(values)
  purrr::map_chr(distinct_values, transform)[match(values, distinct_values)]
}

parse_condition_labels <- function(conditions, cache = NULL) {
  purrr::map_dfr(unique(conditions), function(condition) {
    key <- paste0("label_", digest::digest(condition, algo = "xxhash64"))
    result <- if (is.null(cache)) NULL else cache$get(key)
    if (is.null(result) || inherits(result, "key_missing")) {
      result <- tibble(
        condition = condition,
        receptor_guess_parsed = extract_receptors(condition),
        treatment_guess_parsed = extract_treatment(condition),
        passage_guess_parsed = extract_passage_from_condition(condition)
      )
      if (!is.null(cache)) cache$set(key, result)
    }
    result
  })
}

add_factor_guesses <- function(df, cache = NULL) {
  df %>%
    left_join(parse_condition_labels(df$condition, cache), by = "condition") %>%
    mutate(
      passage_guess_parsed   = if_else(is.na(passage_guess_parsed), passage, passage_guess_parsed),
      
      receptor_guess = coalesce(receptor_pm, receptor_guess_parsed),
      treatment_guess = coalesce(treatment_pm, treatment_guess_parsed),
      passage_guess = coalesce(passage_pm, passage_guess_parsed),
      
      receptor_guess = map_distinct_chr(receptor_guess, canonicalize_receptor_one),
      treatment_guess = map_distinct_chr(treatment_guess, canonicalize_treatment_one),
      passage_guess = map_distinct_chr(passage_guess, standardize_passage_label),
      
      guess_key = make_guess_key(receptor_guess, treatment_guess, passage_guess)
    )
}

parse_conditions_hitl <- function(raw_tbl, cache = NULL) {
  raw_tbl %>%
    distinct(condition_id, .keep_all = TRUE) %>%
    add_factor_guesses(cache) %>%
    distinct(condition_id, file, plate_id, well_id, condition, expt_pm, cell_line_pm,
             receptor_guess, treatment_guess, passage_guess, guess_key) %>%
    mutate(
      expt = coalesce(expt_pm, ""),
      cell_line = coalesce(cell_line_pm, ""),
      receptor = receptor_guess,
      treatment = treatment_guess,
      passage = passage_guess,
      original_guess = paste(receptor_guess, treatment_guess, passage_guess, sep = " | ")
    ) %>%
    update_editor_matching_wells() %>%
    arrange(file, well_id, condition)
}

update_editor_matching_wells <- function(tbl) {
  tbl %>%
    mutate(
      well_id = as.character(well_id),
      receptor = clean_user_factor(receptor),
      treatment = clean_user_factor(treatment),
      passage = map_distinct_chr(passage, standardize_passage_label)
    ) %>%
    group_by(plate_id, passage, receptor, treatment) %>%
    mutate(n_matching_wells = n_distinct(well_id[well_id != ""])) %>%
    ungroup()
}

make_zicht_export_df <- function(df) {
  df %>%
    transmute(
      plate_id = as.character(plate_id),
      well = toupper(str_squish(as.character(well_id))),
      hormone = as.character(treatment),
      receptor = as.character(receptor),
      passage = as.character(passage),
      expt = as.character(expt),
      cell_line = as.character(cell_line)
    ) %>%
    filter(!is.na(well), well != "") %>%
    distinct() %>%
    arrange(plate_id, well)
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
      well      = extract_well_id(well),
      hormone   = as.character(hormone),
      receptor  = as.character(receptor),
      passage   = if ("passage" %in% names(.)) as.character(passage) else NA_character_,
      expt      = if ("expt" %in% names(.)) as.character(expt) else NA_character_,
      cell_line = if ("cell_line" %in% names(.)) as.character(cell_line) else NA_character_,
      receptor_pm  = canonicalize_receptor_combo(receptor),
      treatment_pm = canonicalize_treatment_combo(hormone),
      passage_pm   = if_else(
        is.na(passage) | str_squish(passage) == "",
        NA_character_,
        purrr::map_chr(passage, standardize_passage_label)
      )
    ) %>%
    select(any_of("plate_id"), well, cell_line, expt, receptor_pm, treatment_pm, passage_pm) %>%
    { if (any(.$well == "")) stop("Plate map wells must be A1 through H12."); . }
}

infer_plate_size <- function(wells) {
  wells <- unique(na.omit(wells[wells != ""]))
  if (length(wells) == 0) return(NA_character_)
  
  rows <- str_extract(wells, "^[A-Z]")
  cols <- suppressWarnings(as.integer(str_extract(wells, "[0-9]+$")))
  
  n_rows <- max(match(rows, LETTERS), na.rm = TRUE)
  n_cols <- max(cols, na.rm = TRUE)
  n_total <- n_rows * n_cols
  
  case_when(
    n_rows <= 2 && n_cols <= 3  ~ "6-well",
    n_rows <= 3 && n_cols <= 4  ~ "12-well",
    n_rows <= 4 && n_cols <= 6  ~ "24-well",
    n_rows <= 6 && n_cols <= 8  ~ "48-well",
    n_rows <= 8 && n_cols <= 12 ~ "96-well",
    TRUE ~ paste0(n_total, "-well (inferred)")
  )
}

make_plate_preview_tables <- function(df) {
  df <- df %>%
    mutate(
      well_id  = toupper(str_squish(as.character(well_id))),
      plate_id = as.character(plate_id),
      label = if_else(
        well_id == "" | is.na(well_id),
        NA_character_,
        paste0(
          "<div style='line-height:1.35; padding:4px; background-color:",
          if_else(n_matching_wells > 1, "#FFB81C", "transparent"),
          "'>",
          "<div><strong>Well:</strong> ", well_id, "</div>",
          "<div><strong>Receptor:</strong> ", as.character(receptor), "</div>",
          "<div><strong>Treatment:</strong> ", as.character(treatment), "</div>",
          "<div><strong>Passage:</strong> ", as.character(passage), "</div>",
          "<div><strong>Matching wells:</strong> ", n_matching_wells, "</div>",
          "</div>"
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
      
      list(plate_id = plate_name, table = grid)
    })
}

# ---------------------------
# Math helpers
# ---------------------------
mask_spikes_neighbor_mad <- function(df,
                                     value_col = "value_norm",
                                     threshold = 8) {
  
  y <- df[[value_col]]
  
  df$spike_flag <- FALSE
  df$value_masked <- y
  df$robust_z <- NA_real_
  
  n <- length(y)
  
  if (n < 5)
    return(df)
  
  expected <- rep(NA_real_, n)
  
  expected[2:(n - 1)] <-
    (y[1:(n - 2)] + y[3:n]) / 2
  
  resid <- abs(y - expected)
  
  med <- median(resid, na.rm = TRUE)
  mad_val <- mad(resid, na.rm = TRUE)
  
  if (is.na(mad_val) || mad_val == 0)
    return(df)
  
  robust_z <- (resid - med) / (1.4826 * mad_val)
  
  spike_idx <- which(abs(robust_z) > threshold)
  
  # never automatically remove first/last point
  spike_idx <- setdiff(spike_idx, c(1, n))
  
  if (length(spike_idx) > 0) {
    df$spike_flag[spike_idx] <- TRUE
    df$value_masked[spike_idx] <- NA_real_
  }
  
  df$robust_z <- robust_z
  
  df
}

ols_control_adjust <- function(df, sig_col, ctl_col, fallback = "missing") {
  x <- df[[ctl_col]]
  y <- df[[sig_col]]
  ok <- is.finite(x) & is.finite(y)
  
  df$ols_invalid <- sum(ok) < 3 || length(unique(x[ok])) < 2
  if (df$ols_invalid[1]) {
    df$value_norm <- if (fallback == "raw") ifelse(is.finite(y), y, NA_real_) else NA_real_
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

auc_trapz <- function(x, y, gap_policy = "missing") {
  if (length(x) < 2 || any(!is.finite(x))) return(NA_real_)
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  if (any(diff(x) <= 0)) return(NA_real_)
  if (gap_policy == "missing" && any(!is.finite(y))) return(NA_real_)
  if (gap_policy == "bridge") {
    valid <- is.finite(y)
    x <- x[valid]
    y <- y[valid]
    if (length(x) < 2) return(NA_real_)
  }
  valid_intervals <- is.finite(head(y, -1)) & is.finite(tail(y, -1))
  if (!any(valid_intervals)) return(NA_real_)
  sum((diff(x) * (head(y, -1) + tail(y, -1)) / 2)[valid_intervals])
}

normalize_trajectory <- function(df, signal, control, method = "ratio",
                                 zero_policy = "missing", control_floor = 1,
                                 baseline = FALSE, baseline_policy = "missing",
                                 ols_fallback = "missing") {
  df <- arrange(df, elapsed)
  signal_values <- df[[signal]]
  control_values <- if (control %in% names(df)) df[[control]] else rep(NA_real_, nrow(df))
  df$zero_control <- is.finite(control_values) & control_values == 0
  df$invalid_ratio <- FALSE
  df$ols_invalid <- FALSE
  df$baseline_invalid <- FALSE
  df$trajectory_excluded <- FALSE
  if (method == "none") {
    df$value_norm <- signal_values
  } else if (method == "ols_adj") {
    df <- ols_control_adjust(df, signal, control, ols_fallback)
  } else {
    adjusted_control <- control_values
    if (zero_policy == "floor") {
      if (!is.finite(control_floor) || control_floor <= 0) stop("Control floor must be a positive finite number.")
      adjusted_control[df$zero_control] <- control_floor
    }
    values <- signal_values / adjusted_control
    if (method == "log2ratio") values <- suppressWarnings(log2(values))
    df$invalid_ratio <- !is.finite(values)
    df$value_norm <- ifelse(is.finite(values), values, NA_real_)
    if (zero_policy == "exclude" && any(df$zero_control | df$invalid_ratio)) {
      df$value_norm <- NA_real_
      df$trajectory_excluded <- TRUE
    }
  }
  df$value_norm[!is.finite(df$value_norm)] <- NA_real_
  if (baseline) {
    baseline_value <- df$value_norm[1]
    df$baseline_invalid <- !is.finite(baseline_value) || baseline_value == 0
    if (df$baseline_invalid[1] && baseline_policy == "first_valid") {
      candidates <- which(is.finite(df$value_norm) & df$value_norm != 0)
      baseline_value <- if (length(candidates)) df$value_norm[candidates[1]] else NA_real_
    }
    if (is.finite(baseline_value) && baseline_value != 0) {
      df$value_norm <- df$value_norm / baseline_value
    } else if (baseline_policy != "raw") {
      df$value_norm <- NA_real_
    }
  }
  df$value_norm[!is.finite(df$value_norm)] <- NA_real_
  df
}

# ---------------------------
# Plot/export helpers
# ---------------------------
preview_tabular_file <- function(path, n = 20) {
  lines <- readLines(path, warn = FALSE)
  header_i <- find_data_header_row(lines)
  
  txt <- if (is.na(header_i)) {
    paste(lines, collapse = "\n")
  } else {
    paste(lines[header_i:length(lines)], collapse = "\n")
  }
  
  read_delim(
    I(txt),
    delim = "\t",
    show_col_types = FALSE,
    n_max = n,
    col_types = cols(.default = col_character())
  )
}

clipboard_csv_text <- function(df) {
  if (is.null(df) || nrow(df) == 0) return("")
  paste(capture.output(write.csv(df, row.names = FALSE, na = "")), collapse = "\n")
}

clipboard_matrix_csv_text <- function(x) {
  if (is.null(x) || length(x) == 0) return("")
  paste(capture.output(write.table(x, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE, na = "")), collapse = "\n")
}

classify_treatment_group <- function(treatment) {
  trt <- toupper(as.character(treatment))
  ligands <- c()
  
  if (str_detect(trt, "\\bE2\\b")) ligands <- c(ligands, "E2")
  if (str_detect(trt, "\\bP4\\b")) ligands <- c(ligands, "P4")
  if (str_detect(trt, "\\bDHT\\b")) ligands <- c(ligands, "DHT")
  if (str_detect(trt, "\\b4-OHT\\b")) ligands <- c(ligands, "4-OHT")
  if (str_detect(trt, "\\bRU-486\\b")) ligands <- c(ligands, "RU-486")
  if (str_detect(trt, "\\bFUL\\b")) ligands <- c(ligands, "FUL")
  if (str_detect(trt, "\\bDEX\\b|\\bCORT\\b|\\bGLUCO\\b")) ligands <- c(ligands, "Glucocorticoid")
  if (str_detect(trt, "\\bNET\\b")) ligands <- c(ligands, "NET")
  
  if (length(ligands) == 0) return(as.character(treatment))
  paste(ligands, collapse = " + ")
}

treatment_levels_master <- c(
  "VEH", "E2", "P4", "DHT", "NET", "4-OHT", "RU-486", "FUL", "Glucocorticoid",
  "E2 + P4", "E2 + DHT", "E2 + NET", "E2 + 4-OHT", "E2 + RU-486",
  "E2 + FUL",
  "P4 + DHT", "P4 + 4-OHT", "P4 + RU-486",
  "DHT + 4-OHT", "DHT + RU-486", "4-OHT + RU-486"
)

treatment_color_values <- c(
  "VEH" = "#000000",
  "E2" = "#FB0280",
  "P4" = "#FD8008",
  "DHT" = "#0F80FF",
  "NET" = "#00C896",
  "4-OHT" = "#7A3CFF",
  "RU-486" = "#9E9E9E",
  "FUL" = "#A04E2A",
  "Glucocorticoid" = "#00A878",
  "E2 + P4" = "#C23B8E",
  "E2 + DHT" = "#8A4DFF",
  "E2 + NET" = "#5BB8A0",
  "E2 + 4-OHT" = "#B04DFF",
  "E2 + RU-486" = "#B85A8A",
  "E2 + FUL" = "#C66B52",
  "P4 + DHT" = "#7F9CFF",
  "P4 + 4-OHT" = "#C06A88",
  "P4 + RU-486" = "#C08A60",
  "DHT + 4-OHT" = "#4F5BFF",
  "DHT + RU-486" = "#6B8BB8",
  "4-OHT + RU-486" = "#8E6BA8"
)

compute_auc_export_dims <- function(df) {
  n_receptors <- dplyr::n_distinct(df$receptor)
  n_treatments <- dplyr::n_distinct(df$treatment)
  
  width_mm <- max(89, min(70 + 12 * n_receptors + 6 * n_treatments, 240))
  height_mm <- max(70, 30 + 4 * n_treatments)
  
  list(width_mm = width_mm, height_mm = height_mm)
}

make_treatment_styles <- function(treatments) {
  labels <- sort(unique(as.character(treatments)))
  labels <- labels[!is.na(labels)]
  if (!length(labels)) return(list(linetype = character(), shape = numeric(), color = character()))
  patterns <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash",
                as.vector(outer(1:15, 1:15, function(mark, space) sprintf("%X%X", mark, space))))
  symbols <- c(16, 17, 15, 18, 3, 4, 0:2, 5:14, 19:25)
  colors <- grDevices::hcl.colors(length(labels), "Dark 3")
  groups <- purrr::map_chr(labels, classify_treatment_group)
  unique_groups <- !duplicated(groups) & !duplicated(groups, fromLast = TRUE)
  established <- unique_groups & groups %in% names(treatment_color_values)
  colors[established] <- unname(treatment_color_values[groups[established]])
  list(linetype = setNames(rep_len(patterns, length(labels)), labels),
       shape = setNames(rep_len(symbols, length(labels)), labels),
       color = setNames(colors, labels))
}

make_auc_matrix <- function(data) {
  treatments <- sort(unique(as.character(data$treatment)))
  if (!nrow(data)) return(matrix(c("receptor", "replicate"), nrow = 1))
  columns <- tibble(treatment = treatments, column_key = paste0("treatment_", seq_along(treatments)))
  wide <- data %>%
    mutate(treatment = as.character(treatment)) %>%
    left_join(columns, by = "treatment") %>%
    select(receptor, replicate, column_key, auc) %>%
    pivot_wider(names_from = column_key, values_from = auc) %>%
    arrange(receptor, replicate) %>%
    select(receptor, replicate, all_of(columns$column_key))
  rbind(c("receptor", "replicate", treatments), as.matrix(wide))
}

# ---------------------------
# UI
# ---------------------------
ui <- fluidPage(
  tags$head(
    tags$script(HTML(
      "Shiny.addCustomMessageHandler('copy-to-clipboard', async function(message) {
        const text = (message && message.text) || '';
        const fallback = function() {
          const area = document.createElement('textarea');
          area.value = text;
          area.setAttribute('readonly', '');
          area.style.position = 'fixed';
          area.style.opacity = '0';
          document.body.appendChild(area);
          area.focus();
          area.select();
          document.execCommand('copy');
          document.body.removeChild(area);
        };
        try {
          if (navigator.clipboard && window.isSecureContext) {
            await navigator.clipboard.writeText(text);
          } else {
            fallback();
          }
          Shiny.setInputValue('clipboard_copy_status', { ok: true, nonce: Date.now() }, { priority: 'event' });
        } catch (err) {
          try {
            fallback();
            Shiny.setInputValue('clipboard_copy_status', { ok: true, nonce: Date.now() }, { priority: 'event' });
          } catch (fallbackErr) {
            Shiny.setInputValue('clipboard_copy_status', { ok: false, nonce: Date.now() }, { priority: 'event' });
          }
        }
      });"
    ))
  ),
  titlePanel("Incucyte Multi-File Import + Channel Normalization"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput(
        "files",
        "Upload tab-delimited Incucyte exports (.txt/.tsv/.csv)",
        multiple = TRUE,
        accept = c(".txt", ".tsv", ".csv")
      ),
      actionButton("clear_files", "Clear uploaded files", class = "btn-warning"),
      br(), br(),
      
      checkboxInput("drop_stderr", "Drop '(Std Err ...)' columns", value = TRUE),
      uiOutput("channel_map_ui"),
      
      hr(),
      uiOutput("norm_ui"),
      
      actionButton("run", "Import + Process", class = "btn-primary"),
      
      hr(),
      h4("Downloads"),
      downloadButton("download_prism", "AUC matrix (csv)"),
      downloadButton("download_auc_details", "AUC values and coverage (csv)"),
      downloadButton("download_timecourse", "Timecourse data (csv)"),
      br(), br(),
      actionButton("copy_auc", "Copy AUC to clipboard"),
      actionButton("copy_timecourse", "Copy timelapse to clipboard")
    ),
    
    mainPanel(
      tabsetPanel(
        id = "main_tabs",

        tabPanel(
          "Import",
          h4("Import and normalization checks"),
          uiOutput("numeric_warnings"),
          tags$p("Check plate assignments and channels before Import + Process. Review numerical warnings here after processing."),
          h4("Imported file preview"),
          DTOutput("preview_files_dt")
        ),
        
        tabPanel(
          "Plate map",
          br(),
          fileInput(
            "platemap",
            "Upload plate map (.csv/.zicht)",
            multiple = FALSE,
            accept = c(".csv", ".zicht")
          ),
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
          "Factor editor",
          br(),
          tags$p(
            tags$strong("How it works: "),
            "Each row is a single imported condition. Assign experiment (expt), cell line, receptor, treatment, and passage consistently across a well's channel files, then Apply edits. Plate and well identities are retained internally."
          ),
          fluidRow(
            column(4, actionButton("apply_editor", "Apply edits", class = "btn-primary")),
            column(4, actionButton("reset_editor", "Reset all edits")),
            column(4, downloadButton("download_zicht", "Export .zicht"))
            
          ),
          br(),
          tags$p("Double-click a cell to edit it, or use the controls above the table to update selected rows or all rows. Update all rows includes rows on other pages and rows hidden by filters. Edits affect analysis only after Apply edits."),
          fluidRow(
            column(4, selectInput("bulk_field", "Field", c("Experiment" = "expt", "Cell line" = "cell_line",
                                                            "Passage" = "passage", "Receptor" = "receptor", "Treatment" = "treatment"))),
            column(4, textInput("bulk_value", "Value")),
            column(4,
                   actionButton("bulk_edit", "Update selected rows"),
                   actionButton("bulk_edit_all", "Update all rows"))
          ),
          br(),
          DTOutput("editor_table")
        ),
        
        tabPanel(
          "Plot",
          uiOutput("metadata_warnings"),
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
          ),
          h4("AUC coverage"),
          DTOutput("auc_coverage")
        )
      )
    )
  )
)

# ---------------------------
# Server
# ---------------------------
server <- function(input, output, session) {
  session_cache <- cachem::cache_mem(max_size = 64 * 1024^2)
  label_cache <- cachem::cache_mem(max_size = 8 * 1024^2)
  parsed_files <- new.env(parent = emptyenv())
  session$onSessionEnded(function() {
    session_cache$reset()
    label_cache$reset()
  })
  
  plot_tab_active <- reactive({
    identical(input$main_tabs, "Plot")
  })
  
  uploaded_files_rv <- reactiveVal(NULL)
  preview_files_rv <- reactiveVal(NULL)
  upload_sequence <- reactiveVal(0L)
  plate_choices_rv <- reactiveVal(character())
  
  observeEvent(input$files, {
    req(input$files)
    
    combined <- uploaded_files_rv()
    for (file_index in seq_len(nrow(input$files))) {
      file <- input$files$name[file_index]
      path <- input$files$datapath[file_index]
      parsed <- tryCatch(read_incucyte_file(path, drop_stderr = FALSE), error = function(error) {
        showNotification(paste(file, conditionMessage(error)), type = "error", duration = NULL)
        NULL
      })
      if (is.null(parsed)) next
      upload_sequence(upload_sequence() + 1L)
      file_id <- paste0("file_", upload_sequence())
      parsed_files[[file_id]] <- parsed
      source_key <- plate_source_key(file, list(meta = parsed$meta, dat = select(parsed$dat, -matches("Std Err"))))
      channel_default <- guess_channel(file)
      matches <- if (is.null(combined)) tibble(plate_id = character()) else filter(combined, auto_key == source_key)
      available <- if (nrow(matches)) {
        setdiff(unique(matches$plate_id), matches$plate_id[matches$channel_default == channel_default])
      } else character()
      plate_id <- if (length(available)) available[1] else paste0(
        "plate_", digest::digest(paste(source_key, n_distinct(matches$plate_id) + 1L),
                                algo = "sha256", serialize = FALSE)
      )
      occurrence <- if (is.null(combined)) 1L else sum(combined$file == file) + 1L
      separate_id <- paste0("plate_", digest::digest(paste(source_key, file, occurrence, "separate"),
                                                    algo = "sha256", serialize = FALSE))
      combined <- bind_rows(combined, tibble(
        file = file, path = path, file_id = file_id,
        auto_key = source_key, channel_default = channel_default,
        plate_id = plate_id, separate_id = separate_id
      ))
    }
    if (is.null(combined)) return()
    defaults <- combined %>% distinct(plate_id, .keep_all = TRUE)
    choices <- c(
      setNames(defaults$plate_id, paste0("Plate ", seq_len(nrow(defaults)), " — ", defaults$file)),
      setNames(combined$separate_id, paste0("Separate plate — ", combined$file_id, " — ", combined$file))
    )
    plate_choices_rv(choices)
    uploaded_files_rv(combined)
  })
  
  observeEvent(input$clear_files, {
    uploaded_files_rv(NULL)
    preview_files_rv(NULL)
    plate_choices_rv(character())
    rm(list = ls(parsed_files), envir = parsed_files)
    session_cache$reset()
    label_cache$reset()
    editor_rv(NULL)
    editor_buffer_rv(NULL)
    applied_editor_rv(NULL)
    editor_initialized(FALSE)
  })
  
  observeEvent(uploaded_files_rv(), {
    files_df <- uploaded_files_rv()
    if (is.null(files_df) || nrow(files_df) == 0) {
      preview_files_rv(NULL)
      return()
    }
    
    previews <- purrr::imap(files_df$path, function(path, i) {
      dat <- head(parsed_files[[files_df$file_id[i]]]$dat, 20)
      dat %>% mutate(`..file` = files_df$file[i], `..upload` = files_df$file_id[i], .before = 1)
    })
    
    preview_files_rv(bind_rows(previews))
  }, ignoreNULL = FALSE)
  
  output$channel_map_ui <- renderUI({
    req(uploaded_files_rv())
    files_df <- uploaded_files_rv()
    fns <- files_df$file
    
    tagList(
      h4("Assign channels and plates"),
      tags$p(tags$small("Files are retained until cleared or the app reloads.")),
      tags$p(tags$small("Matching channels should use the same plate. Check automatic suggestions; use Separate plate for independent experiments.")),
      lapply(seq_along(fns), function(i) {
        fname <- fns[i]
        channel_input <- paste0("chan_", files_df$file_id[i])
        plate_input <- paste0("plate_", files_df$file_id[i])
        default <- isolate(input[[channel_input]]) %||% files_df$channel_default[i]
        selected_plate <- isolate(input[[plate_input]]) %||% files_df$plate_id[i]
        
        fluidRow(
          column(8, tags$small(fname)),
          column(4, selectInput(channel_input, NULL, c("GFP", "NIR", "Orange", "Red", "Other"), selected = default)),
          column(12, selectInput(plate_input, "Plate assignment", choices = plate_choices_rv(), selected = selected_plate))
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
      file_id = files_df$file_id,
      plate_id = purrr::map_chr(seq_len(nrow(files_df)), function(i) {
        input[[paste0("plate_", files_df$file_id[i])]] %||% files_df$plate_id[i]
      }),
      channel = purrr::map_chr(seq_len(nrow(files_df)), function(i) {
        val <- input[[paste0("chan_", files_df$file_id[i])]]
        
        if (is.null(val) || is.na(val) || val == "") {
          fname <- files_df$file[i]
          if (str_detect(tolower(fname), "red|nir")) return("NIR")
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
    preview <- plate_map_tbl()
    if ("plate_id" %in% names(preview)) {
      preview <- preview %>% mutate(plate = paste("Plate", match(plate_id, unique(plate_id)))) %>%
        select(plate, everything(), -plate_id)
    }
    preview
  }, striped = TRUE)
  
  output$preview_files_dt <- renderDT({
    req(preview_files_rv())
    datatable(preview_files_rv(), options = list(pageLength = 20, scrollX = TRUE), rownames = FALSE)
  }, server = TRUE)
  
  output$norm_ui <- renderUI({
    req(channel_map())
    
    chans <- sort(unique(channel_map()$channel))
    if (length(chans) == 0) return(NULL)
    
    tagList(
      h4("Normalization"),
      selectInput("signal_channel", "Signal channel", choices = chans, selected = chans[1]),
      selectInput(
        "control_channel",
        "Control channel",
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
          "OLS-adjusted (control-adjusted)" = "ols_adj"
        ),
        selected = "ratio"
      ),
      checkboxInput("baseline_norm", "Baseline-normalize each well trajectory", value = FALSE),
      selectInput("zero_policy", "Zero / invalid control handling",
                  c("Mark undefined points missing" = "missing",
                    "Exclude affected well trajectory" = "exclude",
                    "Replace zero controls with a chosen floor" = "floor")),
      conditionalPanel("input.zero_policy == 'floor'",
                       numericInput("control_floor", "Positive floor (control units)", value = 1, min = 1e-12),
                       tags$p("Choose a floor justified by your measurement scale; this changes ratios.")),
      selectInput("baseline_policy", "Invalid initial baseline",
                  c("Mark trajectory missing" = "missing",
                    "Use first finite nonzero baseline" = "first_valid",
                    "Keep trajectory without baseline scaling" = "raw")),
      selectInput("ols_fallback", "Constant / insufficient OLS control",
                  c("Mark trajectory missing" = "missing", "Keep raw signal" = "raw")),
      selectInput("auc_gap_policy", "AUC with missing points",
                  c("Mark AUC missing" = "missing", "Integrate adjacent valid intervals only" = "segments",
                    "Bridge missing points explicitly" = "bridge")),
      checkboxInput("mask_spikes", "Mask local trajectory spikes", value = FALSE),
      
      conditionalPanel(
        condition = "input.mask_spikes == true",
        numericInput(
          "spike_z_threshold",
          "Spike threshold: robust z-score",
          value = 8,
          min = 3,
          max = 100,
          step = 1
        )
      )
    )
  })
  
  imported_data <- eventReactive(list(input$run, input$clear_files), {
    cm <- channel_map()
    
    imported <- purrr::pmap_dfr(cm, function(file, path, file_id, plate_id, channel) {
      result <- parsed_files[[file_id]]
      req(result)
      meta <- result$meta
      dat_raw <- result$dat
      if (isTRUE(input$drop_stderr)) dat_raw <- select(dat_raw, -matches("Std Err"))
      
      default_passage <- if (!is.na(meta$passage[[1]]) && meta$passage[[1]] != "") {
        standardize_passage_label(meta$passage[[1]])
      } else {
        "Passage_NA"
      }
      
      long <- dat_raw %>%
        pivot_longer(
          cols = -c(datetime, elapsed),
          names_to = "condition",
          values_to = "value"
        ) %>%
        mutate(
          condition = as.character(condition),
          value = suppressWarnings(as.numeric(value)),
          elapsed = suppressWarnings(as.numeric(elapsed)),
          datetime = as.character(datetime),
          file = file,
          plate_id = plate_id,
          channel = channel,
          well_id = toupper(str_squish(extract_well_id(condition))),
          well_key = if_else(well_id != "", well_id, condition),
          replicate_id = paste(plate_id, well_key, sep = " || "),
          passage = default_passage,
          vessel_name = meta$vessel_name[[1]] %||% NA_character_,
          metric = meta$metric[[1]] %||% NA_character_,
          cell_type = meta$cell_type[[1]] %||% NA_character_,
          analysis = meta$analysis[[1]] %||% NA_character_,
          condition_id = paste(file_id, condition, sep = " || "),
          receptor_pm = NA_character_,
          treatment_pm = NA_character_,
          passage_pm = NA_character_,
          cell_line_pm = NA_character_,
          expt_pm = NA_character_
        )
      
      if (!is.null(input$platemap)) {
        long <- join_plate_map(long, plate_map_tbl(), unique(cm$plate_id))
      }
      
      long %>%
        select(
          condition_id, file, plate_id, well_id, well_key, replicate_id, channel, passage,
          vessel_name, metric, cell_type, analysis,
          cell_line_pm, expt_pm,
          receptor_pm, treatment_pm, passage_pm,
          datetime, elapsed, condition, value
        )
    })
    validate(need(all(is.finite(imported$elapsed)), "Some elapsed times are not numeric. Correct the source time column before processing."))
    duplicates <- imported %>% count(plate_id, well_key, channel, elapsed) %>% filter(n > 1)
    validate(need(nrow(duplicates) == 0,
                  "Duplicate plate/well/channel/time observations: assign independent files to separate plates. No wells have been averaged."))
    imported
  }, ignoreInit = TRUE)

  raw_long_auto_base <- reactive({
    req(uploaded_files_rv())
    imported_data()
  })
  
  hitl_default <- reactive({
    req(raw_long_auto_base())
    parse_conditions_hitl(raw_long_auto_base(), label_cache)
  })
  
  editor_rv <- reactiveVal(NULL)
  editor_buffer_rv <- reactiveVal(NULL)
  applied_editor_rv <- reactiveVal(NULL)
  editor_initialized <- reactiveVal(FALSE)
  editor_proxy <- dataTableProxy("editor_table", session = session)
  editable_fields <- c("expt", "cell_line", "passage", "receptor", "treatment")

  editor_display <- function(data) {
    data %>% select(condition_id, file, well_id, expt, cell_line, passage, receptor, treatment,
                    condition, n_matching_wells, original_guess)
  }
  
  observeEvent(hitl_default(), {
    editor_rv(hitl_default())
    editor_buffer_rv(hitl_default())
    applied_editor_rv(NULL)
    editor_initialized(TRUE)
  })
  
  observeEvent(input$reset_editor, {
    req(hitl_default())
    editor_rv(hitl_default())
    editor_buffer_rv(hitl_default())
    applied_editor_rv(NULL)
  })
  
  output$editor_table <- renderDT({
    req(editor_initialized())
    data <- isolate(editor_display(editor_buffer_rv()))
    datatable(data, rownames = FALSE, selection = "multiple",
              editable = list(target = "cell", disable = list(columns = which(!names(data) %in% editable_fields) - 1L)),
              options = list(pageLength = 20, scrollX = TRUE, stateSave = FALSE,
                             columnDefs = list(list(targets = 0, visible = FALSE, searchable = FALSE))))
  }, server = TRUE)

  observeEvent(editor_buffer_rv(), {
    req(editor_initialized())
    replaceData(editor_proxy, editor_display(editor_buffer_rv()), rownames = FALSE,
                resetPaging = FALSE, clearSelection = "none")
  }, ignoreNULL = TRUE)

  update_editor_cells <- function(rows, field, value) {
    data <- editor_buffer_rv()
    if (is.null(data) || !field %in% editable_fields || length(value) != 1L) return()
    rows <- unique(as.integer(rows))
    rows <- rows[!is.na(rows) & rows >= 1L & rows <= nrow(data)]
    if (!length(rows)) return()
    value <- if (field == "passage") standardize_passage_label(value) else clean_user_factor(value)
    data[[field]][rows] <- value
    editor_buffer_rv(data)
  }

  observeEvent(input$editor_table_cell_edit, {
    edit <- input$editor_table_cell_edit
    req(editor_buffer_rv())
    column <- as.integer(edit$col) + 1L
    columns <- names(editor_display(editor_buffer_rv()))
    if (length(column) != 1L || is.na(column) || column < 1L || column > length(columns)) return()
    update_editor_cells(edit$row, columns[column], edit$value)
  })

  observeEvent(input$bulk_edit, {
    req(input$editor_table_rows_selected, input$bulk_field)
    update_editor_cells(input$editor_table_rows_selected, input$bulk_field, input$bulk_value %||% "")
  })

  observeEvent(input$bulk_edit_all, {
    req(editor_buffer_rv(), input$bulk_field)
    update_editor_cells(seq_len(nrow(editor_buffer_rv())), input$bulk_field, input$bulk_value %||% "")
  })
  
  observeEvent(input$apply_editor, {
    req(editor_buffer_rv())
    
    df <- editor_buffer_rv() %>%
      mutate(
        receptor = clean_user_factor(receptor),
        treatment = clean_user_factor(treatment),
        passage = purrr::map_chr(passage, standardize_passage_label),
        factor_key = make_factor_key(receptor, treatment)
      ) %>%
      update_editor_matching_wells()
    
    applied_editor_rv(df)
    editor_rv(df)
    editor_buffer_rv(df)
  }, ignoreInit = TRUE)
  
  current_editor_map <- reactive({
    req(editor_rv())
    
    df <- if (!is.null(applied_editor_rv())) applied_editor_rv() else editor_rv()
    
    df %>%
      mutate(
        receptor = clean_user_factor(receptor),
        treatment = clean_user_factor(treatment),
        passage = purrr::map_chr(passage, standardize_passage_label),
        factor_key = make_factor_key(receptor, treatment)
      ) %>%
      select(condition_id, passage, receptor, treatment, factor_key, expt, cell_line)
  })
  
  zicht_export_df <- reactive({
    req(raw_long_auto_base(), current_editor_map())
    
    condition_metadata() %>%
      distinct(plate_id, well_id, receptor, treatment, passage, expt, cell_line) %>%
      make_zicht_export_df() %>%
      { validate(need(!anyDuplicated(.[c("plate_id", "well")]),
                      "Resolve conflicting metadata between channel files before exporting the plate map.")); . }
  })
  
  output$plate_check_summary <- renderPrint({
    req(raw_long_auto_base())
    
    wells_in_data <- raw_long_auto_base() %>% pull(well_id) %>% unique()
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
      
      if (length(unmatched_data) > 0) cat("Data-only wells:", paste(sort(unmatched_data), collapse = ", "), "\n")
      if (length(unmatched_pm) > 0) cat("Plate-map-only wells:", paste(sort(unmatched_pm), collapse = ", "), "\n")
    }
  })
  
  condition_metadata <- reactive({
    req(raw_long_auto_base(), current_editor_map())
    
    parsed_raw <- raw_long_auto_base() %>% distinct(condition_id, .keep_all = TRUE) %>%
      select(condition_id, plate_id, well_id, well_key, replicate_id, file, condition)
    
    em <- current_editor_map() %>%
      select(condition_id, passage, receptor, treatment, factor_key, expt, cell_line)
    
    parsed_raw %>%
      left_join(em, by = "condition_id", suffix = c("_raw", "")) %>%
      mutate(
        passage = factor(clean_user_factor(passage)),
        receptor = factor(clean_user_factor(receptor)),
        treatment = factor(clean_user_factor(treatment)),
        expt = coalesce(clean_user_factor(expt), ""),
        cell_line = coalesce(clean_user_factor(cell_line), ""),
        factor_key = make_factor_key(receptor, treatment)
      )
  })

  well_metadata <- reactive({
    metadata <- condition_metadata() %>%
      distinct(plate_id, well_key, passage, receptor, treatment, factor_key, expt, cell_line)
    conflicts <- metadata %>% count(plate_id, well_key) %>% filter(n > 1)
    validate(need(nrow(conflicts) == 0,
                  "Channel files disagree on well metadata. Match experiment, cell line, passage, receptor and treatment in the factor editor, then Apply edits."))
    metadata
  })

  receptor_choices <- reactiveVal(character())
  treatment_choices <- reactiveVal(character())
  treatment_styles <- reactive(make_treatment_styles(treatment_choices()))
  observeEvent(condition_metadata(), {
    receptor_choices(sort(unique(as.character(condition_metadata()$receptor))))
    treatment_choices(sort(unique(as.character(condition_metadata()$treatment))))
  })

  plate_layout_data <- reactive({
    req(identical(input$main_tabs, "Plate map"))
    condition_metadata() %>%
      distinct(plate_id, well_id, passage, receptor, treatment) %>%
      group_by(plate_id, passage, receptor, treatment) %>%
      mutate(n_matching_wells = n_distinct(well_id[well_id != ""])) %>%
      ungroup() %>%
      make_plate_preview_tables()
  })
  output$plate_check_layout <- renderUI({
    plate_tables <- plate_layout_data()
    if (length(plate_tables) == 0) return(tags$p("No well-based layout available."))
    tagList(
      purrr::imap(plate_tables, function(plate, plate_name) {
        data <- plate$table
        tagList(
          tags$h4(paste("Plate", match(plate_name, unique(condition_metadata()$plate_id)))),
          tags$table(class = "table table-striped table-bordered",
                     tags$thead(tags$tr(lapply(names(data), tags$th))),
                     tags$tbody(lapply(seq_len(nrow(data)), function(row) {
                       tags$tr(lapply(data[row, ], function(cell) tags$td(HTML(as.character(cell)))))
                     }))),
          tags$br()
        )
      })
    )
  })
  
  wide_joined_passage <- reactive({
    req(raw_long_auto_base())
    raw_long_auto_base() %>%
      select(plate_id, well_id, well_key, replicate_id, channel, elapsed, value) %>%
      pivot_wider(names_from = channel, values_from = value)
  })

  measurement_key <- reactive(digest::digest(wide_joined_passage(), algo = "xxhash64"))

  normalization_settings <- reactive(list(
    method = input$norm_method %||% "ratio", signal = input$signal_channel, control = input$control_channel,
    zero_policy = input$zero_policy %||% "missing", control_floor = input$control_floor %||% 1,
    baseline = isTRUE(input$baseline_norm), baseline_policy = input$baseline_policy %||% "missing",
    ols_fallback = input$ols_fallback %||% "missing", mask_spikes = isTRUE(input$mask_spikes),
    spike_threshold = input$spike_z_threshold %||% 3
  ))
  cached_normalization <- reactive({
    req(wide_joined_passage())
    
    w <- wide_joined_passage()
    method <- input$norm_method %||% "ratio"
    sig <- input$signal_channel
    ctl <- input$control_channel
    
    req(sig, ctl)
    validate(need(sig %in% names(w), "Reprocess after changing channel assignments."),
             need(method == "none" || ctl %in% names(w), "Selected control channel is absent. Check assignments and reprocess."))
    out <- w %>%
      group_by(plate_id, well_key) %>%
      group_modify(~ normalize_trajectory(
        .x, signal = sig, control = ctl, method = method,
        zero_policy = input$zero_policy %||% "missing", control_floor = input$control_floor %||% 1,
        baseline = isTRUE(input$baseline_norm), baseline_policy = input$baseline_policy %||% "missing",
        ols_fallback = input$ols_fallback %||% "missing"
      )) %>%
      ungroup()
    
    if (isTRUE(input$mask_spikes)) {
      threshold <- input$spike_z_threshold %||% 3
      
      out <- out %>%
        group_by(plate_id, well_key) %>%
        arrange(elapsed, .by_group = TRUE) %>%
        group_modify(~ mask_spikes_neighbor_mad(
          .x,
          value_col = "value_norm",
          threshold = threshold
        )) %>%
        ungroup() %>%
        mutate(value_norm = value_masked) %>%
        select(-value_masked)
    } else {
      out <- out %>% mutate(spike_flag = FALSE)
    }
    
    out
  }) %>% bindCache(measurement_key(), normalization_settings(), cache = session_cache)

  normalized_passage <- reactive({
    req(uploaded_files_rv())
    cached_normalization()
  })
  
  stats_long <- reactive({
    req(normalized_passage())
    
    normalized_passage() %>%
      select(plate_id, well_id, well_key, replicate_id, elapsed, value_norm, spike_flag) %>%
      left_join(well_metadata(), by = c("plate_id", "well_key")) %>%
      mutate(plate_label = paste("Plate", match(plate_id, unique(raw_long_auto_base()$plate_id)))) %>%
      arrange(plate_id, well_key, elapsed)
  })

  output$metadata_warnings <- renderUI({
    req(condition_metadata())
    metadata <- condition_metadata() %>%
      distinct(plate_id, well_key, well_id, expt, cell_line, passage, receptor, treatment)
    missing <- metadata %>% filter(
      expt == "" | cell_line == "" | is.na(passage) | passage == "Passage_NA" |
        is.na(receptor) | receptor == "" | is.na(treatment) | treatment == "" | well_id == ""
    )
    repeated <- metadata %>% count(expt, cell_line, passage, receptor, treatment) %>% filter(n > 1)
    conflicts <- metadata %>% count(plate_id, well_key) %>% filter(n > 1)
    messages <- character()
    if (nrow(missing)) messages <- c(messages, paste(nrow(missing), "well metadata records are incomplete. Assign experiment, cell line, passage, receptor and treatment; check unrecognized wells."))
    if (nrow(repeated)) messages <- c(messages, paste(nrow(repeated), "metadata combinations identify multiple wells or plates. Confirm intentional replicates or add missing metadata. Each well remains separate."))
    if (nrow(conflicts)) messages <- c(messages, paste(nrow(conflicts), "wells have conflicting metadata across channels. Resolve these in the factor editor before analysis."))
    if (!length(messages)) return(NULL)
    tags$div(class = "alert alert-warning", role = "alert", tags$strong("Metadata needs review"),
             tags$ul(lapply(messages, tags$li)))
  })

  output$numeric_warnings <- renderUI({
    req(raw_long_auto_base())
    raw <- raw_long_auto_base()
    control_rows <- raw %>% filter(channel == input$control_channel)
    zero_count <- sum(is.finite(control_rows$value) & control_rows$value == 0)
    messages <- character()
    if (zero_count) messages <- c(messages, paste(zero_count, "zero control measurements found. Choose missing points, exclude the affected trajectory, or set a justified positive floor."))
    if (sum(!is.finite(raw$value))) messages <- c(messages, paste(sum(!is.finite(raw$value)), "input measurements are missing or nonfinite."))
    if (identical(input$signal_channel, input$control_channel) && input$norm_method != "none") {
      messages <- c(messages, "Signal and control are the same channel. Select different channels or use None to retain raw signal.")
    }
    normalized <- tryCatch(normalized_passage(), error = function(error) NULL)
    if (is.null(normalized)) {
      messages <- c(messages, "Normalization is unavailable. Check channel assignments and matching metadata, then reprocess if assignments changed.")
    } else {
      if (any(normalized$invalid_ratio)) messages <- c(messages, paste(sum(normalized$invalid_ratio), "ratios are undefined (missing control, division by zero, or invalid log ratio)."))
      if (any(normalized$ols_invalid)) messages <- c(messages, paste(n_distinct(normalized$replicate_id[normalized$ols_invalid]), "well trajectories have constant controls or fewer than three finite pairs; the selected OLS fallback is applied."))
      if (any(normalized$baseline_invalid)) messages <- c(messages, paste(n_distinct(normalized$replicate_id[normalized$baseline_invalid]), "well trajectories have zero or invalid initial baselines; the selected baseline option is applied."))
      if (any(normalized$trajectory_excluded)) messages <- c(messages, "Affected trajectories have been excluded from numerical results and retained as missing values for traceability.")
      if (any(!is.finite(normalized$value_norm))) messages <- c(messages, "Missing normalized points remain. Review the AUC gap policy and coverage table before export.")
    }
    if (!length(messages)) return(tags$div(class = "alert alert-success", "No zero controls or numerical problems detected with the current settings."))
    tags$div(class = "alert alert-warning", role = "alert", tags$strong("Numerical handling needs review"),
             tags$ul(lapply(messages, tags$li)))
  })
  
  output$plot_time_ui <- renderUI({
    req(plot_tab_active())
    req(normalized_passage())
    
    df <- normalized_passage() %>% filter(is.finite(elapsed))
    if (nrow(df) == 0) return(NULL)
    
    rng <- range(df$elapsed, na.rm = TRUE)
    
    sliderInput(
      "plot_time_range",
      "Elapsed time range",
      min = floor(rng[1]),
      max = ceiling(rng[2]),
      value = pmin(ceiling(rng[2]), pmax(floor(rng[1]), isolate(input$plot_time_range) %||% c(floor(rng[1]), ceiling(rng[2])))),
      step = 1
    )
  })
  
  output$plot_receptor_ui <- renderUI({
    req(plot_tab_active())
    levs <- receptor_choices()
    req(length(levs))
    
    checkboxGroupInput("plot_receptors", "Receptors to show", choices = levs,
                       selected = intersect(levs, isolate(input$plot_receptors) %||% levs))
  })
  
  output$plot_treatment_ui <- renderUI({
    req(plot_tab_active())
    levs <- treatment_choices()
    req(length(levs))
    
    checkboxGroupInput("plot_treatments", "Treatments to show", choices = levs,
                       selected = intersect(levs, isolate(input$plot_treatments) %||% levs))
  })
  
  filtered_stats_long <- reactive({
    req(stats_long())
    
    df <- stats_long()
    
    if (!is.null(input$plot_time_range)) {
      df <- df %>%
        filter(elapsed >= input$plot_time_range[1], elapsed <= input$plot_time_range[2])
    }
    
    if (!is.null(input$plot_receptors)) {
      df <- df %>% filter(as.character(receptor) %in% input$plot_receptors)
    }
    
    if (!is.null(input$plot_treatments)) {
      df <- df %>% filter(as.character(treatment) %in% input$plot_treatments)
    }
    
    df
  })
  
  timecourse_export_long <- reactive({
    req(filtered_stats_long())
    
    filtered_stats_long() %>%
      mutate(
        hormone = factor(as.character(treatment), ordered = TRUE),
        receptor = factor(as.character(receptor)),
        passage = factor(str_replace(as.character(passage), "^Passage_", "p")),
        elapsed_hour = elapsed,
        id = replicate_id
      ) %>%
      group_by(plate_id, well_key, cell_line, hormone, receptor, expt, passage, elapsed_hour, id) %>%
      summarise(
        mean_count = mean(value_norm, na.rm = TRUE),
        median_count = median(value_norm, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(cell_line, hormone, receptor, expt, passage, elapsed_hour)
  })
  
  timecourse_prism_export_matrix <- reactive({
    req(filtered_stats_long())
    
    export_long <- filtered_stats_long() %>%
      mutate(
        condition_label = paste(receptor, treatment, passage, expt, cell_line, plate_label, well_key, sep = " | ")
      ) %>%
      arrange(condition_label, elapsed, factor_key)
    
    replicate_map <- export_long %>%
      distinct(condition_label, factor_key) %>%
      group_by(condition_label) %>%
      arrange(factor_key, .by_group = TRUE) %>%
      mutate(rep_idx = row_number()) %>%
      ungroup()
    
    export_wide <- export_long %>%
      left_join(replicate_map, by = c("condition_label", "factor_key")) %>%
      mutate(col_key = paste0(condition_label, "__rep", rep_idx)) %>%
      select(elapsed, col_key, value_norm) %>%
      distinct() %>%
      pivot_wider(names_from = col_key, values_from = value_norm) %>%
      arrange(elapsed)
    
    col_template <- replicate_map %>%
      count(condition_label, name = "n_rep") %>%
      group_by(condition_label) %>%
      summarise(max_rep = max(n_rep), .groups = "drop") %>%
      mutate(col_keys = purrr::map2(condition_label, max_rep, ~ paste0(.x, "__rep", seq_len(.y)))) %>%
      pull(col_keys) %>%
      unlist()
    
    export_wide <- export_wide %>%
      select(elapsed, any_of(col_template))
    
    value_cols <- names(export_wide)[-1]
    condition_header <- c("", str_replace(value_cols, "__rep\\d+$", ""))
    rep_header <- c("elapsed", str_extract(value_cols, "rep\\d+$"))
    
    rbind(condition_header, rep_header, as.matrix(export_wide))
  })
  
  auc_values <- reactive({
    data <- normalized_passage()
    if (!is.null(input$plot_time_range)) {
      data <- data %>% filter(elapsed >= input$plot_time_range[1], elapsed <= input$plot_time_range[2])
    }
    data %>%
      arrange(elapsed) %>%
      group_by(plate_id, well_key) %>%
      summarize(
        auc = auc_trapz(elapsed, value_norm, input$auc_gap_policy %||% "missing"),
        valid_points = sum(is.finite(value_norm)),
        total_points = n(),
        valid_intervals = sum(is.finite(head(value_norm, -1)) & is.finite(tail(value_norm, -1))),
        total_intervals = max(n() - 1L, 0L),
        start_time = min(elapsed), end_time = max(elapsed),
        gap_policy = input$auc_gap_policy %||% "missing",
        .groups = "drop"
      )
  }) %>% bindCache(measurement_key(), normalization_settings(), input$plot_time_range,
                   input$auc_gap_policy %||% "missing", cache = session_cache)

  auc_preview <- reactive({
    filtered_stats_long() %>%
      distinct(plate_id, plate_label, well_key, replicate_id, expt, cell_line, receptor, treatment, passage) %>%
      inner_join(auc_values(), by = c("plate_id", "well_key")) %>%
      arrange(plate_id, plate_label, well_key, replicate_id, expt, cell_line, receptor, treatment, passage)
  })

  output$auc_coverage <- renderDT({
    req(plot_tab_active())
    data <- auc_preview() %>% select(plate_label, well_key, expt, cell_line, passage, receptor, treatment,
                            auc, valid_points, total_points, valid_intervals, total_intervals,
                            start_time, end_time, gap_policy)
    datatable(data, rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  }, server = TRUE)
  
  auc_export_rows <- reactive({
    auc_preview() %>%
      arrange(receptor, treatment, expt, cell_line, passage, plate_id, well_key) %>%
      group_by(receptor, treatment) %>%
      mutate(replicate = row_number()) %>%
      ungroup()
  })

  prism_auc_export_matrix <- reactive({
    make_auc_matrix(auc_export_rows())
  })

  auc_details_export <- reactive({
    auc_export_rows() %>% select(-plate_id, -plate_label, -well_key, -replicate_id)
  })
  
  output$plot <- renderPlot({
    req(filtered_stats_long())
    
    df <- filtered_stats_long()
    styles <- treatment_styles()
    
    if (all(is.na(df$value_norm)) || nrow(df) == 0) {
      plot.new()
      text(0.5, 0.5, "No data remain after current filters.")
      return()
    }
    
    df_ok <- df
    df_spike <- df %>% filter(spike_flag)
    
    p <- ggplot(df_ok, aes(x = elapsed, y = value_norm, group = replicate_id, color = receptor,
                          linetype = treatment, shape = treatment)) +
      geom_line(alpha = 0.6, na.rm = TRUE) +
      geom_point(size = 1, alpha = 0.65, na.rm = TRUE) +
      scale_linetype_manual(values = styles$linetype) +
      scale_shape_manual(values = styles$shape) +
      facet_wrap(~ plate_label + passage) +
      labs(
        x = "Elapsed time",
        y = "Value (normalized)",
        color = "Receptor", linetype = "Treatment", shape = "Treatment",
        title = "Normalized trajectories"
      ) +
      theme_minimal()
    
    if (nrow(df_spike) > 0) {
      p <- p +
        geom_point(
          data = df_spike,
          shape = 4,
          color = "red",
          size = 2,
          aes(group = factor_key)
        )
    }
    
    p
  }) %>% bindCache({ req(plot_tab_active()); filtered_stats_long() }, treatment_styles(), cache = session_cache)
  
  auc_plot_obj <- reactive({
    req(auc_preview())
    
    df <- auc_preview()
    styles <- treatment_styles()
    
    if (nrow(df) == 0 || all(!is.finite(df$auc))) return(NULL)
    
    ggplot(df, aes(x = receptor, y = auc, color = treatment, shape = treatment, group = treatment)) +
      geom_point(position = position_dodge(width = 0.8), size = 2.8, alpha = 0.7, stroke = 0.6, na.rm = TRUE) +
      scale_color_manual(values = styles$color, drop = TRUE) +
      scale_shape_manual(values = styles$shape, drop = TRUE) +
      labs(
        x = "Receptor",
        y = "AUC",
        color = "Treatment", shape = "Treatment",
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
  }) %>% bindCache(auc_preview(), treatment_styles(), cache = session_cache)
  
  output$auc_plot <- renderPlot({
    
    p <- auc_plot_obj()
    
    if (is.null(p)) {
      plot.new()
      text(0.5, 0.5, "No AUC values for current filters.")
      return()
    }
    
    p
  }) %>% bindCache({ req(plot_tab_active()); auc_preview() }, treatment_styles(), cache = session_cache)
  
  
  output$download_prism <- downloadHandler(
    filename = function() paste0("incucyte_prism_auc_", Sys.Date(), ".csv"),
    content = function(file) {
      write.table(
        prism_auc_export_matrix(),
        file = file,
        sep = ",",
        row.names = FALSE,
        col.names = FALSE,
        quote = TRUE,
        na = ""
      )
    }
  )
  
  output$download_auc_details <- downloadHandler(
    filename = function() paste0("incucyte_auc_coverage_", Sys.Date(), ".csv"),
    content = function(file) {
      write.csv(auc_details_export(), file, row.names = FALSE, na = "")
    }
  )

  output$download_zicht <- downloadHandler(
    filename = function() {
      base_name <- uploaded_files_rv()$file[1] %||% "incucyte_platemap"
      paste0(tools::file_path_sans_ext(base_name), ".zicht")
    },
    content = function(file) {
      write.table(
        zicht_export_df(),
        file = file,
        sep = ",",
        row.names = FALSE,
        col.names = TRUE,
        quote = TRUE,
        na = ""
      )
    }
  )
  
  output$download_timecourse <- downloadHandler(
    filename = function() paste0("incucyte_timecourse_", Sys.Date(), ".csv"),
    content = function(file) {
      write.table(
        timecourse_prism_export_matrix(),
        file = file,
        sep = ",",
        row.names = FALSE,
        col.names = FALSE,
        quote = TRUE,
        na = ""
      )
    }
  )
  
  observeEvent(input$copy_auc, {
    session$sendCustomMessage(
      "copy-to-clipboard",
      list(text = clipboard_matrix_csv_text(prism_auc_export_matrix()))
    )
  })
  
  observeEvent(input$copy_timecourse, {
    session$sendCustomMessage(
      "copy-to-clipboard",
      list(text = clipboard_matrix_csv_text(timecourse_prism_export_matrix()))
    )
  })
  
  observeEvent(input$clipboard_copy_status, {
    if (isTRUE(input$clipboard_copy_status$ok)) {
      showNotification("Copied to clipboard.", type = "message", duration = 2)
    } else {
      showNotification("Clipboard copy failed in this browser session.", type = "error", duration = 4)
    }
  })
  
  output$download_auc_plot_png <- downloadHandler(
    filename = function() paste0("incucyte_auc_plot_", Sys.Date(), ".png"),
    content = function(file) {
      p <- isolate({
        if (isTRUE(plot_tab_active())) auc_plot_obj() else NULL
      })
      
      if (is.null(p)) {
        png(file, width = 1800, height = 1200, res = 300)
        plot.new()
        text(0.5, 0.5, "Open the Plot tab first to render the AUC plot.")
        dev.off()
        return()
      }
      
      df <- auc_preview()
      
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
      p <- isolate({
        if (isTRUE(plot_tab_active())) auc_plot_obj() else NULL
      })
      
      if (is.null(p)) {
        svg(filename = file, width = 3.5, height = 2.75, bg = "white")
        plot.new()
        text(0.5, 0.5, "Open the Plot tab first to render the AUC plot.")
        dev.off()
        return()
      }
      
      df <- auc_preview()
      
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
}

shinyApp(ui, server)
