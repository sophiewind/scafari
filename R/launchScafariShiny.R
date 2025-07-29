#' Launch Scafari Shiny Application
#'
#' Launches the Scafari Shiny application for data visualization.
#'
#' @examples
#' if (interactive()) {
#'     launchScafariShiny()
#' }
#'
#' @export
#' @return shiny app
launchScafariShiny <- function() {
    shiny::shinyApp(ui = app_ui(), server = app_server)
}
