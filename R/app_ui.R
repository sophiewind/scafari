app_ui <- function() {
    navbarPage(
        useShinyjs(), # Enable shinyjs
        createUploadUI(),  # Upload tab
        createSequencingUI(),  # Seq tab
        
        # Panel tab ------------------------------------------------------------
        tabPanel("Panel",
            fluidPage(
                conditionalPanel(condition = "output.file_ready",
                    createPanelUI(),),
                conditionalPanel(condition = "!output.file_ready",
                    fluidRow(
                        style = "text-align: center; margin-top: 50px;",
                        p(style = "color: black;",
                            "Please upload a file to view sequencing
                            information."))), hr(),)),

        # Variants tab ---------------------------------------------------------
        tabPanel("Variants",
            fluidPage(wellPanel(fluidRow(
                    createFilteringUI())),
                conditionalPanel(
                    condition = "output.plots_visible == true",
                    createVariantUI(),))),

        # Explore Variants tab -------------------------------------------------
        tabPanel("Explore profiles",
            fluidPage(createExploreVariantSelectionUI(),
                conditionalPanel(
                    condition = "output.continue == true",
                    h1("Identify cell clusters"),
                    h2("Select number of clusters"),
                    withLoader(plotOutput("kneeplot")), loader = "dnaspin",
                    fluidRow(
                        numericInput("n_clust", "Number of clusters:",
                            value = 3, min = 0, max = 10, step = 1),
                        div(style = "display:inline-block; float:center",
                            actionButton("kmeans_btn",
                                label = "Start kmeans",
                                icon = icon("play"),
                                style = paste0("color: #fff; ",
                                    "background-color: #337ab7; ",
                                    "border-color: #2e6da4"))),
                        hr()),
                    createExploreVariantUI())),
            hr(),
        )
    )
}
