createExploreVariantSelectionUI <- function() {
  tagList(
    fluidRow(
      sidebarLayout(
        sidebarPanel(width = 3,
                     tags$style(HTML(".dataTable {font-size: 12px; overflow-y: scroll; }")),
                     dataTableOutput("data_table_var2"),  
                     hr(),
                     actionButton('submit_var', label = 'Select variants', icon = icon('play-circle'), 
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),),
        
        mainPanel(h1('Explore profiles'),
                  plotOutput('hm_1', height = "800px" ),
                  
        )
      ),
      div(style = "display:inline-block; float:center", 
          actionButton('continue_var', 'Continue with this variant selection', icon = icon('thumbs-up'), 
                       style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
      hr()
    )
  )
}