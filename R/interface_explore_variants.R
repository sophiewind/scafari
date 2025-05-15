createExploreVariantUI <- function() {
  tags$head(
    tags$style(HTML("
    .shiny-output-error-validation {
      color: red;
    }
  "))
  )
  
  tagList(
    useShinyjs(),
    conditionalPanel(
      condition = "output.plots_visible_2 == true",
      
      h2('Cluster cells by VAF'),
      h4('Variants included in Clustering:'),
      uiOutput("selected_rows_2"),
      h4('Clustering method used:'),
      uiOutput("selected_method"),
      
      fluidRow(
        withLoader(
          plotOutput("cluster_plot"),  
          loader = 'dnaspin'
        )
      ),
      
     conditionalPanel(
        condition = "output.plots_visible_3 == true",
      
        h2('\n'),
        fluidRow(withLoader(plotOutput("vaf_hm", height = "800px"), 
                            loader = 'dnaspin')),
        hr(),
        
        h2('Explore variant profiles'),
        h3('VAF in clusters'),
        fluidRow(withLoader(plotOutput("vaf_violin"), loader = 'dnaspin')),
        
        h3('VAF distribution map'),
        fluidRow(withLoader(plotOutput("vaf_map"), loader = 'dnaspin')),
        
        h3('Numerical genotype in clusters'),
        fluidRow(withLoader(plotOutput('ana_bar'), loader = 'dnaspin'))
      )
    )
  ) 
}