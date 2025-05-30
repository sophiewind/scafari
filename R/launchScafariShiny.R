#' Launch the Scafari Shiny App
#' 
#' @examples
#' \dontrun{
#' launchScafariShiny()
#' }
#' 
#' @export
#' 
#' @return shiny app
launchScafariShiny <- function(){
  shiny::shinyApp(ui = app_ui(), server = app_server)
}