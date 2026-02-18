# app.R
library(shiny)
library(tidyverse)
library(readr)

# ---- helpers ----
guess_delim <- function(path) {
  # Incucyte exports are often tab-delimited .txt; fall back to csv
  ext <- tools::file_ext(path)
  if (tolower(ext) %in% c("txt", "tsv")) "\t" else ","
}

read_incucyte <- function(path) {
  delim <- guess_delim(path)
  df <- read_delim(path, delim = delim, show_col_types = FALSE)
  
  # --- Heuristic: adapt this to your actual Incucyte export format ---
  # Common patterns:
  # - A "Well" column
  # - A "Time" / "Elapsed" column, or time as the first column
  # - A measurement column (e.g., "Mean Green Object Intensity" or similar)
  #
  # Try to find likely columns:
  nm <- names(df)
  
  well_col <- nm[str_detect(tolower(nm), "^well$|well id|wellid|^wells?$")][1]
  time_col <- nm[str_detect(tolower(nm), "^time$|elapsed|time \\(h\\)|time \\(hr\\)|hours")][1]
  
  # measurement column: first numeric column that's not time
  num_cols <- nm[map_lgl(df, is.numeric)]
  meas_col <- setdiff(num_cols, time_col)[1]
  
  if (is.na(well_col) || is.na(time_col) || is.na(meas_col)) {
    stop("Could not auto-detect well/time/value columns. Update read_incucyte() heuristics.")
  }
  
  df %>%
    transmute(
      well = as.character(.data[[well_col]]),
      time = as.numeric(.data[[time_col]]),
      value = as.numeric(.data[[meas_col]])
    )
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("Incucyte Multi-Channel Import + Normalization"),
  sidebarLayout(
    sidebarPanel(
      fileInput("files", "Upload Incucyte export files", multiple = TRUE,
                accept = c(".txt", ".tsv", ".csv")),
      
      uiOutput("channel_map_ui"),
      hr(),
      
      selectInput("signal_channel", "Signal channel", choices = NULL),
      selectInput("control_channel", "Control channel (for normalization)", choices = NULL),
      
      radioButtons("norm_method", "Normalization method",
                   choices = c("Ratio (signal/control)" = "ratio",
                               "Log2 ratio" = "log2ratio"),
                   selected = "ratio"),
      
      actionButton("run", "Import + Normalize", class = "btn-primary"),
      hr(),
      downloadButton("download_tidy", "Download normalized tidy CSV")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Join checks", verbatimTextOutput("join_checks")),
        tabPanel("Preview (raw long)", tableOutput("preview_raw")),
        tabPanel("Preview (wide)", tableOutput("preview_wide")),
        tabPanel("Preview (normalized)", tableOutput("preview_norm")),
        tabPanel("Plot", plotOutput("plot"))
      )
    )
  )
)

# ---- server ----
server <- function(input, output, session) {
  
  # Build file->channel mapping UI dynamically
  output$channel_map_ui <- renderUI({
    req(input$files)
    file_names <- input$files$name
    
    tagList(
      h4("Assign channels"),
      lapply(seq_along(file_names), function(i) {
        fluidRow(
          column(8, tags$small(file_names[i])),
          column(4,
                 selectInput(
                   inputId = paste0("chan_", i),
                   label = NULL,
                   choices = c("GFP", "NIR", "Orange", "Red", "Other"),
                   selected = ifelse(str_detect(tolower(file_names[i]), "nir"), "NIR",
                                     ifelse(str_detect(tolower(file_names[i]), "gfp|green"), "GFP", "Other"))
                 )
          )
        )
      })
    )
  })
  
  # Collect channel assignments
  channel_map <- reactive({
    req(input$files)
    tibble(
      file = input$files$name,
      path = input$files$datapath,
      channel = map_chr(seq_along(input$files$name), ~ input[[paste0("chan_", .x)]])
    )
  })
  
  # Update signal/control choices after files are present
  observeEvent(channel_map(), {
    chans <- sort(unique(channel_map()$channel))
    updateSelectInput(session, "signal_channel", choices = chans, selected = chans[1])
    updateSelectInput(session, "control_channel", choices = chans, selected = if (length(chans) > 1) chans[2] else chans[1])
  })
  
  # Import files into long tidy format
  raw_long <- eventReactive(input$run, {
    cm <- channel_map()
    
    out <- purrr::pmap_dfr(cm, function(file, path, channel) {
      dat <- read_incucyte(path)
      dat %>%
        mutate(file = file, channel = channel)
    })
    
    out
  })
  
  # Wide for normalization
  wide <- reactive({
    req(raw_long())
    raw_long() %>%
      select(well, time, channel, value) %>%
      group_by(well, time, channel) %>%
      summarize(value = mean(value, na.rm = TRUE), .groups = "drop") %>%  # collapse duplicates
      pivot_wider(names_from = channel, values_from = value)
  })
  
  normalized <- reactive({
    req(wide(), input$signal_channel, input$control_channel)
    w <- wide()
    
    sig <- input$signal_channel
    ctl <- input$control_channel
    
    if (!(sig %in% names(w)) || !(ctl %in% names(w))) {
      stop("Selected signal/control channels not found after pivot. Check channel assignment.")
    }
    
    w %>%
      mutate(
        value_norm = case_when(
          input$norm_method == "ratio" ~ .data[[sig]] / .data[[ctl]],
          input$norm_method == "log2ratio" ~ log2(.data[[sig]] / .data[[ctl]]),
          TRUE ~ NA_real_
        )
      )
  })
  
  # Join checks / QC
  output$join_checks <- renderPrint({
    req(raw_long())
    rl <- raw_long()
    
    counts <- rl %>%
      distinct(channel, well) %>%
      count(channel, name = "n_wells")
    
    times <- rl %>%
      distinct(channel, time) %>%
      count(channel, name = "n_times")
    
    list(
      wells_per_channel = counts,
      times_per_channel = times,
      channels = sort(unique(rl$channel)),
      total_rows = nrow(rl)
    )
  })
  
  # Previews
  output$preview_raw <- renderTable({
    req(raw_long())
    head(raw_long(), 20)
  })
  
  output$preview_wide <- renderTable({
    req(wide())
    head(wide(), 20)
  })
  
  output$preview_norm <- renderTable({
    req(normalized())
    head(normalized(), 20)
  })
  
  # Plot (simple)
  output$plot <- renderPlot({
    req(normalized())
    df <- normalized()
    
    ggplot(df, aes(x = time, y = value_norm, group = well)) +
      geom_line(alpha = 0.3) +
      labs(x = "Time", y = "Normalized value") +
      theme_minimal()
  })
  
  # Download
  output$download_tidy <- downloadHandler(
    filename = function() paste0("incucyte_normalized_", Sys.Date(), ".csv"),
    content = function(file) {
      write_csv(normalized(), file)
    }
  )
}

shinyApp(ui, server)