# app.R
# IncuZicht: Plate Map + Incucyte Import + Join
# Updates implemented:
# - Default plate size is 24-well
# - Plate map upload NEVER crashes the app:
#     * shows error message on invalid file
#     * AUTO-SWITCHES plate size if the uploaded plate map is a different size
# - Keeps prior functionality:
#     * Plate Map editing (long table default, optional grid)
#     * Import up to 3 Incucyte files
#     * Join Incucyte long data to plate map by Well, respecting include == TRUE

library(shiny)
library(rhandsontable)

# -------------------------
# Plate map helpers
# -------------------------
plate_dims <- function(n_wells) {
  switch(
    as.character(n_wells),
    "6"  = list(nrow = 2, ncol = 3),
    "12" = list(nrow = 3, ncol = 4),
    "24" = list(nrow = 4, ncol = 6),
    "48" = list(nrow = 6, ncol = 8),
    "96" = list(nrow = 8, ncol = 12),
    stop("Unsupported plate size: ", n_wells)
  )
}

make_plate_map <- function(n_wells) {
  d <- plate_dims(n_wells)
  rows <- LETTERS[seq_len(d$nrow)]
  cols <- seq_len(d$ncol)
  
  grid <- expand.grid(Row = rows, Col = cols, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$Well <- sprintf("%s%02d", grid$Row, grid$Col)
  
  grid[["include"]] <- TRUE
  
  grid[["hormone a"]] <- ""
  grid[["hormone b"]] <- ""
  grid[["receptor a"]] <- ""
  grid[["receptor b"]] <- ""
  grid[["passage number"]] <- ""
  grid[["experiment number"]] <- ""
  grid[["notes"]] <- ""
  
  grid[, c(
    "Well", "Row", "Col",
    "include",
    "hormone a", "hormone b",
    "receptor a", "receptor b",
    "passage number", "experiment number",
    "notes"
  )]
}

infer_plate_size_from_plate_map <- function(df) {
  # Prefer row count if it matches canonical plates; fallback to Row/Col unique counts
  n <- nrow(df)
  if (n %in% c(6, 12, 24, 48, 96)) return(as.integer(n))
  
  if (!all(c("Row", "Col") %in% names(df))) return(NA_integer_)
  nrow_u <- length(unique(df$Row))
  ncol_u <- length(unique(df$Col))
  
  if (nrow_u == 2 && ncol_u == 3) return(6L)
  if (nrow_u == 3 && ncol_u == 4) return(12L)
  if (nrow_u == 4 && ncol_u == 6) return(24L)
  if (nrow_u == 6 && ncol_u == 8) return(48L)
  if (nrow_u == 8 && ncol_u == 12) return(96L)
  
  NA_integer_
}

validate_plate_map <- function(df, n_wells) {
  req_cols <- c(
    "Well", "Row", "Col", "include",
    "hormone a", "hormone b",
    "receptor a", "receptor b",
    "passage number", "experiment number"
  )
  missing <- setdiff(req_cols, names(df))
  if (length(missing)) stop("Plate map missing required columns: ", paste(missing, collapse = ", "))
  
  d <- plate_dims(n_wells)
  expected <- d$nrow * d$ncol
  if (nrow(df) != expected) stop("Plate map has ", nrow(df), " rows; expected ", expected, " for a ", n_wells, "-well plate.")
  
  if (any(!grepl("^[A-Z][0-9]{2}$", df$Well))) stop("Some 'Well' values are not like A01, B12, etc.")
  
  rows <- LETTERS[seq_len(d$nrow)]
  if (any(!df$Row %in% rows)) stop("Some 'Row' values are not valid for this plate size.")
  if (any(!df$Col %in% seq_len(d$ncol))) stop("Some 'Col' values are not valid for this plate size.")
  
  if (!is.logical(df$include)) {
    df$include <- tolower(as.character(df$include)) %in% c("true", "t", "1", "yes", "y")
  }
  
  if (!("notes" %in% names(df))) df[["notes"]] <- ""
  
  df
}

factor_fields <- c(
  "hormone a", "hormone b",
  "receptor a", "receptor b",
  "passage number", "experiment number"
)

# -------------------------
# Incucyte import helpers
# -------------------------
read_incucyte_txt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  vessel <- NA_character_
  vessel_idx <- grep("^\\s*Vessel Name\\s*:", lines, ignore.case = TRUE)
  if (length(vessel_idx) >= 1) {
    vessel <- sub("^\\s*Vessel Name\\s*:\\s*", "", lines[vessel_idx[1]], ignore.case = TRUE)
    vessel <- trimws(vessel)
    if (!nzchar(vessel)) vessel <- NA_character_
  }
  
  hdr_idx <- grep("^\\s*Date Time\\b", lines, ignore.case = FALSE)
  if (length(hdr_idx) == 0) stop("Could not find a header line starting with 'Date Time' in: ", basename(path))
  
  dat <- read.table(
    text = paste(lines[hdr_idx[1]:length(lines)], collapse = "\n"),
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  
  if (!("Date Time" %in% names(dat))) stop("Missing 'Date Time' column in: ", basename(path))
  if (!("Elapsed" %in% names(dat))) stop("Missing 'Elapsed' column in: ", basename(path))
  
  list(vessel = vessel, wide = dat)
}

normalize_well <- function(x) {
  x <- toupper(trimws(as.character(x)))
  m <- regexec("^([A-Z]+)\\s*0*([0-9]+)$", x)
  reg <- regmatches(x, m)
  out <- x
  ok <- lengths(reg) == 3
  if (any(ok)) {
    row <- vapply(reg[ok], `[[`, character(1), 2)
    col <- as.integer(vapply(reg[ok], `[[`, character(1), 3))
    out[ok] <- sprintf("%s%02d", row, col)
  }
  out
}

incucyte_wide_to_long <- function(wide_df) {
  id_cols <- c("Date Time", "Elapsed")
  well_cols <- setdiff(names(wide_df), id_cols)
  
  norm_cols <- normalize_well(well_cols)
  names(wide_df)[match(well_cols, names(wide_df))] <- norm_cols
  well_cols <- norm_cols
  
  mat <- as.matrix(wide_df[well_cols])
  value <- as.vector(mat)
  
  data.frame(
    `Date Time` = rep(wide_df[["Date Time"]], times = length(well_cols)),
    Elapsed = rep(wide_df[["Elapsed"]], times = length(well_cols)),
    Well = rep(well_cols, each = nrow(wide_df)),
    Value = value,
    stringsAsFactors = FALSE
  )
}

guess_plate_size_from_wells <- function(well_names_A01) {
  rows <- sub("^([A-Z]).*$", "\\1", well_names_A01)
  cols <- suppressWarnings(as.integer(sub("^[A-Z]+", "", well_names_A01)))
  rows <- rows[!is.na(rows) & nzchar(rows)]
  cols <- cols[!is.na(cols)]
  if (!length(rows) || !length(cols)) return(NA_integer_)
  
  nrow_u <- length(unique(rows))
  ncol_u <- length(unique(cols))
  
  if (nrow_u == 2 && ncol_u == 3) return(6L)
  if (nrow_u == 3 && ncol_u == 4) return(12L)
  if (nrow_u == 4 && ncol_u == 6) return(24L)
  if (nrow_u == 6 && ncol_u == 8) return(48L)
  if (nrow_u == 8 && ncol_u == 12) return(96L)
  NA_integer_
}

# -------------------------
# UI
# -------------------------
ui <- fluidPage(
  titlePanel("IncuZicht"),
  tabsetPanel(
    id = "main_tabs",
    
    tabPanel(
      "Plate Map",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            "plate_size", "Plate size",
            choices = c("6" = 6, "12" = 12, "24" = 24, "48" = 48, "96" = 96),
            selected = 24  # DEFAULT 24
          ),
          actionButton("new_map", "New blank plate map"),
          tags$hr(),
          checkboxInput("show_grid", "Show grid view (A–H by 1–12)", value = FALSE),
          selectInput(
            "grid_field",
            "Grid edits write to:",
            choices = c(factor_fields, "notes"),
            selected = "hormone a"
          ),
          tags$hr(),
          fileInput("upload_map", "Upload plate map (.csv)", accept = c(".csv")),
          downloadButton("download_map", "Download plate map (.csv)"),
          tags$hr(),
          helpText("Use include=FALSE to mark empty/unused wells (excluded from later analysis).")
        ),
        mainPanel(
          conditionalPanel(
            condition = "input.show_grid == true",
            h4("Grid view (editable)"),
            p("Grid edits the selected field. To exclude wells, use the long table to toggle include."),
            rHandsontableOutput("plate_grid")
          ),
          conditionalPanel(
            condition = "input.show_grid == false",
            h4("Long table view (editable)"),
            rHandsontableOutput("plate_long")
          ),
          tags$br(),
          verbatimTextOutput("plate_map_status"),
          tags$hr()
        )
      )
    ),
    
    tabPanel(
      "Import Data",
      sidebarLayout(
        sidebarPanel(
          fileInput(
            "incu_files",
            "Upload up to 3 Incucyte data files (.txt)",
            multiple = TRUE,
            accept = c(".txt", ".tsv", ".csv")
          ),
          helpText("Expected: optional 'Vessel Name:' line, then tab-delimited header 'Date Time', 'Elapsed', well columns."),
          tags$hr(),
          checkboxInput("make_long", "Build long format (recommended)", value = TRUE),
          helpText("Wells are normalized to A01 format for joining.")
        ),
        mainPanel(
          h4("Import summary"),
          tableOutput("import_summary"),
          tags$hr(),
          h4("Preview: selected file"),
          selectInput("preview_file", "Choose file to preview:", choices = character(0)),
          h5("Wide (first 8 rows)"),
          tableOutput("preview_wide"),
          conditionalPanel(
            condition = "input.make_long == true",
            h5("Long (first 12 rows)"),
            tableOutput("preview_long")
          )
        )
      )
    ),
    
    tabPanel(
      "Merged Data",
      sidebarLayout(
        sidebarPanel(
          helpText("Incucyte data joined to plate map (by Well)."),
          tags$hr(),
          selectInput("merged_file", "Choose imported file:", choices = character(0)),
          checkboxInput("drop_excluded", "Drop excluded wells (include==FALSE)", value = TRUE),
          checkboxInput("only_mapped", "Keep only wells present in plate map", value = TRUE)
        ),
        mainPanel(
          h4("Join checks"),
          verbatimTextOutput("join_checks"),
          tags$hr(),
          h4("Merged preview (first 20 rows)"),
          tableOutput("merged_preview")
        )
      )
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # Prevents "plate_size change" observer from overwriting a freshly uploaded plate map
  suppress_plate_reset <- reactiveVal(FALSE)
  
  # DEFAULT plate map is 24-well
  plate_map <- reactiveVal(make_plate_map(24))
  
  observeEvent(input$new_map, {
    plate_map(make_plate_map(as.integer(input$plate_size)))
  })
  
  observeEvent(input$plate_size, {
    if (isTRUE(suppress_plate_reset())) {
      suppress_plate_reset(FALSE)
      return(NULL)
    }
    plate_map(make_plate_map(as.integer(input$plate_size)))
  }, ignoreInit = TRUE)
  
  # ---- Plate map upload: safe + auto-switch ----
  observeEvent(input$upload_map, {
    req(input$upload_map)
    
    tryCatch({
      df <- read.csv(input$upload_map$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      
      inferred <- infer_plate_size_from_plate_map(df)
      if (is.na(inferred)) {
        stop("Could not infer plate size from uploaded plate map (row count or Row/Col layout not recognized).")
      }
      
      # Auto-switch if mismatch
      if (as.integer(input$plate_size) != inferred) {
        suppress_plate_reset(TRUE)
        updateSelectInput(session, "plate_size", selected = inferred)
        showNotification(
          paste0("Auto-switched plate size to ", inferred, "-well to match uploaded plate map."),
          type = "warning", duration = 6
        )
      }
      
      df <- validate_plate_map(df, inferred)  # validate against inferred size
      plate_map(df)
      
      showNotification("Plate map uploaded successfully.", type = "message", duration = 4)
      
    }, error = function(e) {
      showNotification(
        paste0("Plate map upload failed: ", conditionMessage(e)),
        type = "error",
        duration = 10
      )
    })
  })
  
  # ---- Long table ----
  output$plate_long <- renderRHandsontable({
    rhandsontable(
      plate_map(),
      rowHeaders = NULL,
      stretchH = "all",
      height = 580
    ) %>%
      hot_col("Well", readOnly = TRUE) %>%
      hot_col("Row", readOnly = TRUE) %>%
      hot_col("Col", readOnly = TRUE) %>%
      hot_col("include", type = "checkbox")
  })
  
  observeEvent(input$plate_long, {
    df <- hot_to_r(input$plate_long)
    if (!is.logical(df$include)) {
      df$include <- tolower(as.character(df$include)) %in% c("true", "t", "1", "yes", "y")
    }
    plate_map(df)
  })
  
  # ---- Grid view ----
  output$plate_grid <- renderRHandsontable({
    df <- plate_map()
    d  <- plate_dims(as.integer(input$plate_size))
    rows <- LETTERS[seq_len(d$nrow)]
    cols <- sprintf("%02d", seq_len(d$ncol))
    
    field <- input$grid_field
    
    mat <- matrix("", nrow = d$nrow, ncol = d$ncol, dimnames = list(rows, cols))
    mat[cbind(match(df$Row, rows), df$Col)] <- as.character(df[[field]])
    
    display_mat <- mat
    excluded <- !df$include
    if (any(excluded)) {
      for (i in which(excluded)) {
        rr <- match(df$Row[i], rows)
        cc <- df$Col[i]
        display_mat[rr, cc] <- paste0("(x) ", display_mat[rr, cc])
      }
    }
    
    grid_df <- as.data.frame(display_mat, stringsAsFactors = FALSE)
    grid_df <- cbind(Row = rownames(grid_df), grid_df)
    rownames(grid_df) <- NULL
    
    rhandsontable(grid_df, stretchH = "all", height = 580) %>%
      hot_col("Row", readOnly = TRUE)
  })
  
  observeEvent(input$plate_grid, {
    req(input$plate_grid)
    grid_df <- hot_to_r(input$plate_grid)
    
    d <- plate_dims(as.integer(input$plate_size))
    rows <- LETTERS[seq_len(d$nrow)]
    cols <- sprintf("%02d", seq_len(d$ncol))
    field <- input$grid_field
    
    df <- plate_map()
    
    for (r in rows) {
      row_i <- which(grid_df$Row == r)
      for (c in seq_len(d$ncol)) {
        col_nm <- cols[c]
        val <- as.character(grid_df[row_i, col_nm, drop = TRUE])
        val <- sub("^\\(x\\)\\s*", "", val)
        df[[field]][df$Row == r & df$Col == c] <- ifelse(is.na(val), "", val)
      }
    }
    plate_map(df)
  })
  
  output$download_map <- downloadHandler(
    filename = function() paste0("plate_map_", input$plate_size, "well.csv"),
    content = function(file) write.csv(plate_map(), file, row.names = FALSE, na = "")
  )
  
  output$plate_map_status <- renderPrint({
    df <- plate_map()
    cat("Plate:", input$plate_size, "well\n")
    cat("Included wells:", sum(df$include), " / ", nrow(df), "\n\n")
    
    nonempty <- vapply(
      c(factor_fields, "notes"),
      function(nm) sum(nzchar(trimws(as.character(df[[nm]])))),
      integer(1)
    )
    cat("Non-empty cells by field:\n")
    print(nonempty)
    
    cat("\nPreview (first 10 wells):\n")
    print(utils::head(df, 10))
  })
  
  # ---- Import data state ----
  incu_data <- reactiveVal(list())
  
  observeEvent(input$incu_files, {
    req(input$incu_files)
    
    if (nrow(input$incu_files) > 3) {
      showNotification("Please upload at most 3 files.", type = "error")
      return(NULL)
    }
    
    files <- input$incu_files
    out <- list()
    
    for (i in seq_len(nrow(files))) {
      path <- files$datapath[i]
      fname <- files$name[i]
      
      parsed <- read_incucyte_txt(path)
      wide <- parsed$wide
      
      wells_raw <- setdiff(names(wide), c("Date Time", "Elapsed"))
      wells_norm <- normalize_well(wells_raw)
      names(wide)[match(wells_raw, names(wide))] <- wells_norm
      
      long <- NULL
      if (isTRUE(input$make_long)) {
        long <- incucyte_wide_to_long(wide)
      }
      
      guessed_plate <- guess_plate_size_from_wells(wells_norm)
      
      out[[fname]] <- list(
        filename = fname,
        vessel = parsed$vessel,
        wide = wide,
        long = long,
        wells = wells_norm,
        n_timepoints = nrow(wide),
        n_wells = length(wells_norm),
        guessed_plate = guessed_plate
      )
    }
    
    incu_data(out)
    
    updateSelectInput(session, "preview_file", choices = names(out), selected = names(out)[1])
    updateSelectInput(session, "merged_file",  choices = names(out), selected = names(out)[1])
    
    # Warn if any file's wells don't match selected plate size (best-effort)
    sel <- as.integer(input$plate_size)
    mism <- vapply(out, function(x) !is.na(x$guessed_plate) && x$guessed_plate != sel, logical(1))
    if (any(mism)) {
      bad <- paste(names(out)[mism], collapse = ", ")
      showNotification(
        paste0(
          "Plate size mismatch: selected ", sel,
          "-well, but these file(s) look like different plates: ", bad,
          ". Joining can still proceed; choose 'Keep only wells present in plate map' to filter."
        ),
        type = "warning",
        duration = 12
      )
    }
  })
  
  output$import_summary <- renderTable({
    dat <- incu_data()
    if (!length(dat)) return(NULL)
    
    data.frame(
      file = vapply(dat, `[[`, character(1), "filename"),
      vessel = vapply(dat, function(x) ifelse(is.na(x$vessel), "", x$vessel), character(1)),
      timepoints = vapply(dat, `[[`, integer(1), "n_timepoints"),
      wells = vapply(dat, `[[`, integer(1), "n_wells"),
      guessed_plate = vapply(dat, function(x) ifelse(is.na(x$guessed_plate), "", as.character(x$guessed_plate)), character(1)),
      stringsAsFactors = FALSE
    )
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$preview_wide <- renderTable({
    dat <- incu_data()
    req(length(dat), input$preview_file)
    utils::head(dat[[input$preview_file]]$wide, 8)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$preview_long <- renderTable({
    req(isTRUE(input$make_long))
    dat <- incu_data()
    req(length(dat), input$preview_file)
    long <- dat[[input$preview_file]]$long
    req(!is.null(long))
    utils::head(long, 12)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  # ---- JOIN: Incucyte long + plate map ----
  merged_data <- reactive({
    dat <- incu_data()
    req(length(dat), input$merged_file)
    
    file_obj <- dat[[input$merged_file]]
    req(!is.null(file_obj$long))
    
    pm <- plate_map()
    
    # (1) plate-map-based filtering if requested
    if (isTRUE(input$only_mapped)) {
      file_long <- file_obj$long[file_obj$long$Well %in% pm$Well, , drop = FALSE]
    } else {
      file_long <- file_obj$long
    }
    
    # (2) join by Well
    merged <- merge(
      file_long,
      pm,
      by = "Well",
      all.x = TRUE,
      sort = FALSE
    )
    
    # Respect include flag
    if (isTRUE(input$drop_excluded)) {
      merged <- merged[isTRUE(merged$include), , drop = FALSE]
    }
    
    merged$source_file <- file_obj$filename
    merged$vessel <- ifelse(is.na(file_obj$vessel), "", file_obj$vessel)
    
    merged
  })
  
  output$join_checks <- renderPrint({
    dat <- incu_data()
    if (!length(dat) || !nzchar(input$merged_file)) {
      cat("Upload data files to enable joining.\n")
      return()
    }
    
    file_obj <- dat[[input$merged_file]]
    pm <- plate_map()
    
    sel_plate <- as.integer(input$plate_size)
    guessed <- file_obj$guessed_plate
    
    cat("Selected plate size:", sel_plate, "\n")
    cat("Guessed plate size from file wells:", ifelse(is.na(guessed), "NA", guessed), "\n\n")
    
    wells_in_file <- unique(file_obj$wells)
    wells_in_pm <- pm$Well
    
    cat("File wells (normalized):", length(wells_in_file), "\n")
    cat("Plate map wells:", length(wells_in_pm), "\n")
    
    extra <- setdiff(wells_in_file, wells_in_pm)
    missing <- setdiff(wells_in_pm, wells_in_file)
    
    cat("\nWells in file but not on this plate map:", length(extra), "\n")
    if (length(extra) && length(extra) <= 20) cat(paste(extra, collapse = ", "), "\n")
    if (length(extra) > 20) cat(paste(c(extra[1:20], "..."), collapse = ", "), "\n")
    
    cat("\nWells on plate map but not in file:", length(missing), "\n")
    if (length(missing) && length(missing) <= 20) cat(paste(missing, collapse = ", "), "\n")
    if (length(missing) > 20) cat(paste(c(missing[1:20], "..."), collapse = ", "), "\n")
    
    cat("\nIncluded wells in plate map:", sum(pm$include), " / ", nrow(pm), "\n")
    cat("Join behavior:\n")
    cat(" - Keep only mapped wells:", isTRUE(input$only_mapped), "\n")
    cat(" - Drop excluded wells:", isTRUE(input$drop_excluded), "\n")
  })
  
  output$merged_preview <- renderTable({
    m <- merged_data()
    utils::head(m, 20)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
}

shinyApp(ui, server)