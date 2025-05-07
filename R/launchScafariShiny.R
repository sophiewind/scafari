#' Launch the Scafari Shiny App
#' 
#' @examples
#' launchScafariShiny()
#' 
#' @export
launchScafariShiny <- function(){
  shiny::shinyApp(ui = app_ui(), server = app_server)
}