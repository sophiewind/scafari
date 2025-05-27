createPanelUI <- function(){
  tagList(
    h2("Panel Information"),  # Header for the panel
    fluidRow(
      column(width=6,
             withLoader(plotOutput("panel_plot1")), loader = 'dnaspin'),
      column(width=6,
             withLoader(plotOutput("panel_plot2")), loader = 'dnaspin')
    ),
    fluidRow(
      h2('Amplicon Overview'),
      withLoader(dataTableOutput("data_table_overview"))
    ),
    fluidRow(
      h2('Amplicon Distribution'),
      withLoader(plotOutput("panel_plot3", height = "800px"), loader = 'dnaspin')
    ),
    fluidRow(
      h2('Amplicon Performance'),
      withLoader(plotOutput("panel_plot4"), loader = 'dnaspin')
    ),
    fluidRow(
      h2('Panel uniformity plot',bsButton("pu", label = "", icon = icon("info"), style = "info", size = "extra-small")),
      textOutput("panelUniformityText"),
      bsPopover(
        id = "pu",
        title = "Panel uniformity",
        content = paste0(
          "The plot indicates the normalized read counts per amplicon. Amplicons that meet at least 20% of the average depth of coverage are plotted in blue others in orange."
        ),
        placement = "right",
        trigger = "hover",
        options = list(container = "body")
      ),
      withLoader(plotlyOutput("panel_plot5"), loader = 'dnaspin')
    )
  )
}
