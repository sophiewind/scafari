createSequencingUI <- function() {
  tagList(
    fluidRow(
      h2('Overview'),
      column(width = 6,
             plotOutput("seq_plot1", height = "300px"),
             withLoader(plotOutput("seq_plot2", height = "300px"), loader = 'dnaspin')
      ),
      column(width = 6,
             withLoader(plotOutput("seq_plot3", height = "600px"), loader = 'dnaspin')
      )
    ),
    hr(),
    fluidRow(
      h2('Sample information'),
      withLoader(dataTableOutput("data_table_sample"), loader = 'dnaspin')
    ),
    hr(),
    fluidRow(
      h2('Sequencing information'),
      withLoader(dataTableOutput("data_table_sequencing"), loader = 'dnaspin')
    ),
    hr(),
    fluidRow(
      h2('Mapping'),
      withLoader(dataTableOutput("data_table_mapping"), loader = 'dnaspin')
    )
  )
}

createNumericInputWithPopover <- function(id, label, popoverContent, ...) {
  tagList(
    numericInput(id, label, ...),
    bsPopover(id, title = label, content = popoverContent, trigger = "hover")
  )
}