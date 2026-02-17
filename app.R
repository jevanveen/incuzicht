# app.R
# IncuZicht: Plate Map + Incucyte Data Import (up to 3 files)

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
  if (nrow(df) != expected) stop("Plate map has ", nrow(df), " rows; expected ", expected, ".")
  
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

# Parse a file like:
# Vessel Name: ...
# Date Time\tElapsed\tA1\tA2...
# ...data rows...
read_incucyte_txt <- function(path) {
  lines <- readLines(path, warn = FALSE)
  
  # Vessel Name (optional)
  vessel <- NA_character_
  vessel_idx <- grep("^\\s*Vessel Name\\s*:", lines, ignore.case = TRUE)
  if (length(vessel_idx) >= 1) {
    vessel <- sub("^\\s*Vessel Name\\s*:\\s*", "", lines[vessel_idx[1]], ignore.case = TRUE)
    vessel <- trimws(vessel)
    if (!nzchar(vessel)) vessel <- NA_character_
  }
  
  # Find header line (starts with "Date Time")
  hdr_idx <- grep("^\\s*Date Time\\b", lines, ignore.case = FALSE)
  if (length(hdr_idx) == 0) stop("Could not find a header line starting with 'Date Time' in: ", basename(path))
  
  # Read table from header onward
  dat <- read.table(
    text = paste(lines[hdr_idx[1]:length(lines)], collapse = "\n"),
    header = TRUE,
    sep = "\t",
    quote = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    comment.char = ""
  )
  
  # Minimal validation
  if (!("Date Time" %in% names(dat))) stop("Missing 'Date Time' column in: ", basename(path))
  if (!("Elapsed" %in% names(dat))) stop("Missing 'Elapsed' column in: ", basename(path))
  
  # Keep Date Time as character for now; you can parse later if you want
  # dat[["Date Time"]] <- as.POSIXct(dat[["Date Time"]], format = "%m/%d/%Y %I:%M:%S %p", tz = "America/Los_Angeles")
  
  list(vessel = vessel, wide = dat)
}

incucyte_wide_to_long <- function(wide_df) {
  id_cols <- c("Date Time", "Elapsed")
  well_cols <- setdiff(names(wide_df), id_cols)
  
  mat <- as.matrix(wide_df[well_cols])
  # column-major: all rows for A1, then all rows for A2, etc.
  value <- as.vector(mat)
  
  data.frame(
    `Date Time` = rep(wide_df[["Date Time"]], times = length(well_cols)),
    Elapsed = rep(wide_df[["Elapsed"]], times = length(well_cols)),
    Well = rep(well_cols, each = nrow(wide_df)),
    Value = value,
    stringsAsFactors = FALSE
  )
}

guess_plate_size_from_wells <- function(well_names) {
  # Supports A1..H12 etc. Returns 6/12/24/48/96 or NA if ambiguous.
  # Assumes well names like A1, A2, ... D6 (no leading zeros) as in your files.  [oai_citation:2‡AM_18_013026_greenintensity_data.txt](sediment://file_0000000001c871fd86381789a035c2aa)
  rows <- sub("^([A-Za-z]).*$", "\\1", well_names)
  cols <- suppressWarnings(as.integer(sub("^[A-Za-z]+", "", well_names)))
  rows <- toupper(rows)
  rows <- rows[!is.na(rows) & nzchar(rows)]
  cols <- cols[!is.na(cols)]
  
  if (!length(rows) || !length(cols)) return(NA_integer_)
  
  nrow <- length(unique(rows))
  ncol <- length(unique(cols))
  
  # Map to canonical plates
  if (nrow == 2 && ncol == 3) return(6L)
  if (nrow == 3 && ncol == 4) return(12L)
  if (nrow == 4 && ncol == 6) return(24L)
  if (nrow == 6 && ncol == 8) return(48L)
  if (nrow == 8 && ncol == 12) return(96L)
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
            selected = 96
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
          helpText("Expected: optional 'Vessel Name:' line, then tab-delimited header 'Date Time', 'Elapsed', well columns (A1, A2, ...)."),
          tags$hr(),
          checkboxInput("make_long", "Also build long format (recommended)", value = TRUE),
          checkboxInput("try_guess_plate", "Guess plate size from well names", value = TRUE)
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
    )
  )
)

# -------------------------
# Server
# -------------------------
server <- function(input, output, session) {
  
  # ---- Plate map state ----
  plate_map <- reactiveVal(make_plate_map(96))
  
  observeEvent(input$new_map, {
    plate_map(make_plate_map(as.integer(input$plate_size)))
  })
  
  observeEvent(input$plate_size, {
    plate_map(make_plate_map(as.integer(input$plate_size)))
  }, ignoreInit = TRUE)
  
  observeEvent(input$upload_map, {
    req(input$upload_map)
    df <- read.csv(input$upload_map$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    df <- validate_plate_map(df, as.integer(input$plate_size))
    plate_map(df)
  })
  
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
      wells <- setdiff(names(wide), c("Date Time", "Elapsed"))
      
      long <- NULL
      if (isTRUE(input$make_long)) {
        long <- incucyte_wide_to_long(wide)
      }
      
      guessed_plate <- NA_integer_
      if (isTRUE(input$try_guess_plate)) {
        guessed_plate <- guess_plate_size_from_wells(wells)
      }
      
      out[[fname]] <- list(
        filename = fname,
        vessel = parsed$vessel,
        wide = wide,
        long = long,
        n_timepoints = nrow(wide),
        n_wells = length(wells),
        guessed_plate = guessed_plate
      )
    }
    
    incu_data(out)
    
    updateSelectInput(
      session, "preview_file",
      choices = names(out),
      selected = names(out)[1]
    )
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
    req(length(dat))
    req(input$preview_file)
    wide <- dat[[input$preview_file]]$wide
    utils::head(wide, 8)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
  
  output$preview_long <- renderTable({
    req(isTRUE(input$make_long))
    dat <- incu_data()
    req(length(dat))
    req(input$preview_file)
    long <- dat[[input$preview_file]]$long
    req(!is.null(long))
    utils::head(long, 12)
  }, striped = TRUE, bordered = TRUE, spacing = "s")
}

shinyApp(ui, server)