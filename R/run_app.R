#' Launch the incuzicht Shiny app
#' @export
run_incuzicht <- function(...) {
  app_dir <- system.file("app", package = "incuzicht")
  
  if (app_dir == "") {
    stop("Could not find the app. Try reinstalling `incuzicht`.", call. = FALSE)
  }
  
  shiny::runApp(app_dir, ...)
}
