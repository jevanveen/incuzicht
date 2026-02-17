# app.R
# Minimal Shiny module to let users CREATE / EDIT / UPLOAD / DOWNLOAD a multi-well plate map
# Supports 6, 12, 24, 48, 96-well plates.

library(shiny)
library(rhandsontable)

# ---- helpers ----
plate_dims <- function(n_wells) {
  switch(
    as.character(n_wells),
    "6"  = list(nrow = 2, ncol = 3),   # A-B x 1-3
    "12" = list(nrow = 3, ncol = 4),   # A-C x 1-4
    "24" = list(nrow = 4, ncol = 6),   # A-D x 1-6
    "48" = list(nrow = 6, ncol = 8),   # A-F x 1-8
    "96" = list(nrow = 8, ncol = 12),  # A-H x 1-12
    stop("Unsupported plate size: ", n_wells)
  )
}

make_plate_map <- function(n_wells) {
  d <- plate_dims(n_wells)
  rows <- LETTERS[seq_len(d$nrow)]
  cols <- seq_len(d$ncol)
  
  grid <- expand.grid(Row = rows, Col = cols, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$Well <- sprintf("%s%02d", grid$Row, grid$Col)
  
  # Common, CRAN-friendly default fields; users can overwrite any cells
  grid$Group     <- ""
  grid$Treatment <- ""
  grid$Dose      <- ""
  grid$Replicate <- ""
  grid$Notes     <- ""
  
  # Order columns nicely
  grid <- grid[, c("Well", "Row", "Col", "Group", "Treatment", "Dose", "Replicate", "Notes")]
  grid
}

validate_plate_map <- function(df, n_wells) {
  req_cols <- c("Well", "Row", "Col")
  missing <- setdiff(req_cols, names(df))
  if (length(missing)) stop("Plate map missing required columns: ", paste(missing, collapse = ", "))
  
  d <- plate_dims(n_wells)
  if (nrow(df) != d$nrow * d$ncol) stop("Plate map has ", nrow(df), " rows; expected ", d$nrow * d$ncol, ".")
  
  # Basic Well format check like A01
  if (any(!grepl("^[A-Z][0-9]{2}$", df$Well))) stop("Some 'Well' values are not like A01, B12, etc.")
  df
}

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
      fileInput("upload_map", "Upload plate map (.csv)", accept = c(".csv")),
      downloadButton("download_map", "Download plate map (.csv)"),
      tags$hr(),
      checkboxInput("show_grid", "Show as well grid (A–H by 1–12)", value = TRUE),
      helpText("Tip: Edit cells directly in the table. Required columns: Well, Row, Col.")
    ),
    mainPanel(
      conditionalPanel(
        condition = "input.show_grid == true",
        h4("Grid view (editable)"),
        rHandsontableOutput("plate_grid")
      ),
      conditionalPanel(
        condition = "input.show_grid == false",
        h4("Long table view (editable)"),
        rHandsontableOutput("plate_long")
      ),
      tags$hr(),
      verbatimTextOutput("plate_map_status")
    )
  )
)

server <- function(input, output, session) {
  
  plate_map <- reactiveVal(make_plate_map(96))
  
  observeEvent(input$new_map, {
    plate_map(make_plate_map(as.integer(input$plate_size)))
  })
  
  observeEvent(input$plate_size, {
    # Keep things simple: changing plate size resets the map
    plate_map(make_plate_map(as.integer(input$plate_size)))
  }, ignoreInit = TRUE)
  
  observeEvent(input$upload_map, {
    req(input$upload_map)
    df <- read.csv(input$upload_map$datapath, stringsAsFactors = FALSE, check.names = FALSE)
    df <- validate_plate_map(df, as.integer(input$plate_size))
    plate_map(df)
  })
  
  # Editable LONG table
  output$plate_long <- renderRHandsontable({
    rhandsontable(
      plate_map(),
      rowHeaders = NULL,
      stretchH = "all",
      height = 520
    )
  })
  
  observeEvent(input$plate_long, {
    plate_map(hot_to_r(input$plate_long))
  })
  
  # Editable GRID view: rows A.., columns 01..; stores Group by default (users can switch to other fields later)
  output$plate_grid <- renderRHandsontable({
    df <- plate_map()
    d  <- plate_dims(as.integer(input$plate_size))
    rows <- LETTERS[seq_len(d$nrow)]
    cols <- sprintf("%02d", seq_len(d$ncol))
    
    # Choose what to display in grid (default: Group)
    mat <- matrix("", nrow = d$nrow, ncol = d$ncol, dimnames = list(rows, cols))
    idx <- match(sprintf("%s%02d", df$Row, df$Col), sprintf("%s%02d", df$Row, df$Col)) # stable
    mat[cbind(match(df$Row, rows), df$Col)] <- df$Group
    
    grid_df <- as.data.frame(mat, stringsAsFactors = FALSE)
    grid_df <- cbind(Row = rownames(grid_df), grid_df)
    rownames(grid_df) <- NULL
    
    rhandsontable(grid_df, stretchH = "all", height = 520) %>%
      hot_col("Row", readOnly = TRUE)
  })
  
  observeEvent(input$plate_grid, {
    # Write grid edits back into plate_map$Group
    grid_df <- hot_to_r(input$plate_grid)
    d <- plate_dims(as.integer(input$plate_size))
    rows <- LETTERS[seq_len(d$nrow)]
    cols <- sprintf("%02d", seq_len(d$ncol))
    
    df <- plate_map()
    # grid_df columns: Row, 01, 02, ...
    for (r in rows) {
      row_i <- which(grid_df$Row == r)
      for (c in seq_len(d$ncol)) {
        col_nm <- cols[c]
        val <- as.character(grid_df[row_i, col_nm, drop = TRUE])
        df$Group[df$Row == r & df$Col == c] <- ifelse(is.na(val), "", val)
      }
    }
    plate_map(df)
  })
  
  output$download_map <- downloadHandler(
    filename = function() {
      paste0("plate_map_", input$plate_size, "well.csv")
    },
    content = function(file) {
      write.csv(plate_map(), file, row.names = FALSE, na = "")
    }
  )
  
  output$plate_map_status <- renderPrint({
    df <- plate_map()
    cat("Current plate map:", nrow(df), "wells\n")
    cat("Non-empty Group entries:", sum(nzchar(df$Group)), "\n")
    cat("Preview:\n")
    print(utils::head(df, 10))
  })
}

shinyApp(ui, server)