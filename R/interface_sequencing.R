createSequencingUI <- function() {
    tabPanel(
        "Sequencing",
        fluidPage(
            # Conditional display based on file upload
            conditionalPanel(
                condition = "output.file_ready",
                tagList(
                    fluidRow(
                        h2("Overview"),
                        column(
                            width = 6,
                            plotOutput("seq_plot1", height = "300px"),
                            withLoader(
                                plotOutput(
                                    "seq_plot2",
                                    height = "300px"
                                ),
                                loader = "dnaspin"
                            )
                        ),
                        column(
                            width = 6,
                            withLoader(
                                plotOutput("seq_plot3",
                                    height = "600px"
                                ),
                                loader = "dnaspin"
                            )
                        )
                    ),
                    hr(),
                    fluidRow(
                        h2("Sample information"),
                        withLoader(dataTableOutput("data_table_sample"),
                            loader = "dnaspin"
                        )
                    ),
                    hr(),
    fluidRow(
      column(width = 4,
             fluidRow(
               h2('Sequencing information'),
               withLoader(dataTableOutput("data_table_sequencing"), loader = 'dnaspin')
             ),
             hr(),
             fluidRow(
               h2('Mapping'),
               withLoader(dataTableOutput("data_table_mapping"), loader = 'dnaspin')
               ),
             ),
      column(width = 8,
             withLoader(plotlyOutput("seq_plot4", height = "600px"), loader = 'dnaspin')
             )
      ),
    hr(),
    fluidRow(
      h2('Tapestri pipeline information'),
      withLoader(dataTableOutput("data_table_tapestri"), loader = 'dnaspin')
            # Message displayed when the file is not uploaded
            conditionalPanel(
                condition = "!output.file_ready",
                fluidRow(
                    style = "text-align: center; margin-top: 50px;",
                    p(
                        style = "color: black;",
                        "Please upload a file to view sequencing
                            information."
                    )
                )
            ),
            hr(),
        )
    )
}

createNumericInputWithPopover <- function(id, label, popoverContent, ...) {
    tagList(
        numericInput(id, label, ...),
        bsPopover(id,
            title = label, content = popoverContent,
            trigger = "hover"
        )
    )
}
