#' Run Shiny app for power estimation
#'
#'
#' @export
#'
runShiny <- function(){

  require(shiny)
  require(plotly)
  require(reshape2)

  appDir <- system.file("shinyApp", package = "scPower")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing the package.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
