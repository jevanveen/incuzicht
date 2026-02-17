# app.R
# IncuZicht: Plate Map Builder (6/12/24/48/96 well)
# Updates requested:
# 1) Long table view is the default (show_grid = FALSE)
# 2) Move the dynamic status output ABOVE the horizontal separator

library(shiny)
library(rhandsontable)

# ---- helpers ----
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

# ---- app ----
ui <- fluidPage(
  titlePanel("IncuZicht: Plate Map Input"),
  sidebarLayout(
    sidebarPanel(
      selectInput(
        "plate_size", "Plate size",
        choices = c("6" = 6, "12" = 12, "24" = 24, "48" = 48, "96" = 96),
        selected = 96
      ),
      actionButton("new_map", "New blank plate map"),
      tags$hr(),
      checkboxInput("show_grid", "Show grid view (A–H by 1–12)", value = FALSE), # <-- default long table
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
      # Table outputs (dynamic based on input.show_grid)
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
      
      # Status moved ABOVE the horizontal separator
      tags$br(),
      verbatimTextOutput("plate_map_status"),
      tags$hr()
    )
  )
)

server <- function(input, output, session) {
  
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
    content = function(file) {
      write.csv(plate_map(), file, row.names = FALSE, na = "")
    }
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
}

shinyApp(ui, server)