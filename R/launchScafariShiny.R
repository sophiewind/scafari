#' Launch the Scafari Shiny App
#'
#' @export
launchScafariShiny <- function(){
  shiny::shinyApp(ui = app_ui(), server = app_server)
}