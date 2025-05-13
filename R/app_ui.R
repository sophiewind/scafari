app_ui <- function() {
  navbarPage(
    useShinyjs(), # Enable shinyjs
    
    # Upload tab -----------------------------------------------------------------
    tabPanel(
      "Upload",
      fluidPage(
        width = 3, # Set sidebar width to 3 columns (out of 12)
        h3("Data Upload"),
        fileInput("upload", "Upload HDF5 File", accept = ".h5"),
        verbatimTextOutput("class_output"),
        verbatimTextOutput("file_contents"),
        div(id = "text"),
        tableOutput("files")
      )
    ),
    
    # Sequencing tab -----------------------------------------------------------
    tabPanel(
      "Sequencing",
      fluidPage(
        # Conditional display based on file upload
        conditionalPanel(
          condition = "output.file_ready",
          createSequencingUI(),
        ),
        # Message displayed when the file is not uploaded
        conditionalPanel(
          condition = "!output.file_ready",
          fluidRow(
            style = "text-align: center; margin-top: 50px;",
            p(
              style = "color: black;",
              "Please upload a file to view sequencing information."
            )
          )
        ),
        hr(),
      )
    ),
    
    # Panel tab ----------------------------------------------------------------
    tabPanel(
      "Panel",
      fluidPage(
        conditionalPanel(
          condition = "output.file_ready",
          createPanelUI(),
        ),
        # Message displayed when the file is not uploaded
        conditionalPanel(
          condition = "!output.file_ready",
          fluidRow(
            style = "text-align: center; margin-top: 50px;",
            p(
              style = "color: black;",
              "Please upload a file to view sequencing information."
            )
          )
        ),
        hr(),
      )
    ),
    
    # Variants tab ---------------------------------------------------------------
    tabPanel(
      "Variants",
      fluidPage(
        wellPanel(fluidRow(
          createFilteringUI()
        )),
        conditionalPanel(
          condition = "output.plots_visible == true",
          createVariantUI(),
        )
      )
    ),
    
    # Explore Variants tab -------------------------------------------------------
    tabPanel(
      "Explore profiles",
      fluidPage(
        createExploreVariantSelectionUI(),
        conditionalPanel(
          condition = "output.continue == true",
          h1("Identify cell clusters"),
          h2("Select clustering method"),
          
          radioButtons(
            "radio",
            "Select option",
            choices = list("k-means" = "kmeans", "leiden" = "leiden", "DBSCAN" = "dbscan"),
            selected = "kmeans"
          ),
          
          h2("Setup clustering parameters"),
          
          # Panel for k-means
          conditionalPanel(
            condition = "input.radio == 'kmeans'",
            fluidRow(
              withLoader(plotOutput("kneeplot"), loader = "dnaspin"),
              numericInput("n_clust", "Number of clusters:",
                           value = 3, min = 0, max = 10, step = 1
              ),
              div(
                style = "display:inline-block; float:center",
                actionButton("kmeans_btn",
                             label = "Start kmeans",
                             icon = icon("play"),
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                )
              ),
              hr()
            )
          ),
          
          # Panel for Leiden
          conditionalPanel(
            condition = "input.radio == 'leiden'",
            fluidRow(
              numericInput("resolution", "Resolution parameter:",
                           value = 0.5, min = 0, max = 1, step = 0.01
              ),
              div(
                style = "display:inline-block; float:center",
                actionButton("kmeans_btn",
                             label = "Start Leiden",
                             icon = icon("play"),
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                )
              ),
              hr()
            )
          ),
          
          # Panel for DBSCAN
          conditionalPanel(
            condition = "input.radio == 'dbscan'",
            fluidRow(
              withLoader(plotOutput("edgeplot"), loader = "dnaspin"),
              
              numericInput("minPts", "minPts:",
                           value = 5, min = 1, max = 100, step = 1
              ),
              numericInput("eps", "eps.value:",
                           value = 0.5, min = 0, max = 5, step = 0.1
              ),
              div(
                style = "display:inline-block; float:center",
                actionButton("kmeans_btn",
                             label = "Start DBSCAN",
                             icon = icon("play"),
                             style = "color: #fff; background-color: #337ab7; border-color: #2e6da4"
                )
              ),
              hr()
            )
          ),
          
          createExploreVariantUI()
        )
      ),
      hr(),
    )
  )
}