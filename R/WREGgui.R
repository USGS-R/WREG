#' Function to open the WREG GUI
#' @examples
#' \dontrun{ 
#' WREGgui()
#' }
#' @export
#' @import shiny
#' @import rmarkdown
#' @importFrom DT renderDataTable
#' @importFrom DT dataTableOutput
WREGgui <- function() {
  appDir <- system.file("shiny", package = "WREG")
  if (appDir == "") {
    stop("Could not find GUI directory. Try re-installing `WQWREG`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal",launch.browser=TRUE)
}