#' Launch the Scafari Shiny App
#' 
#' @examples
#' launchScafariShiny()
#' 
#' @export
#' 
#' @return shiny app
launchScafariShiny <- function(){
  shiny::shinyApp(ui = app_ui(), server = app_server)
}