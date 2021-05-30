#' Run EpiScope shiny app
#'
#' @return The shiny app will open
#'
#' @param dev Run the applicaiton in developer mode
#'
#' @examples
#' \dontrun{
#' EpiScope()
#' }
#' @export

EpiScope <- function(dev=FALSE) {
  appDir <- system.file("shiny", package="EpiScope")
  if (appDir == "") {
    stop("Could not find EpiScope. Try re-installing `EpiScope`.",
         call. = FALSE)
  }
  if (dev) {
    options(shiny.autoreload=TRUE)
  }
  shiny::runApp(appDir, display.mode="normal")
}
