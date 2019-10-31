#' Run Shiny app for power estimation
#'
#' @import shiny
#' @import plotly
#' @import ggplot2
#' @import reshape2
#'
#' @export
#'
runShiny <- function(){

  require(shiny)
  require(plotly)
  require(ggplot2)
  require(reshape2)

  appDir <- system.file("shinyApp", package = "powerScPop")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing the package.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
