createVariantUI <- function() {
  tagList(
    h2("Variant Information"),
    fluidRow(
      column(width=6,
             withLoader(plotOutput("var_plot1", height = "300px")), loader = 'dnaspin'),
      column(width=6,
             withLoader(plotOutput("var_plot2", height = "300px")), loader = 'dnaspin')
    ),
    
    fluidRow(
      h2('Overview Filtered Variants'),
      dataTableOutput("data_table_var")  
    ),
    fluidRow(
      h2('Variant Allele Frequency'),
      column(11,
             withLoader(plotOutput("var_plot3", height = "800px"), loader = 'dnaspin')
      ),
      column(1,
             plotOutput("legend")),
    ),
    fluidRow(
      h2('Genotype of Filtered Variants'),
      column(width=6,
             withLoader(plotOutput("var_plot4")), loader = 'dnaspin'),
      column(width=6,
             withLoader(plotOutput("var_plot5")), loader = 'dnaspin')
    )
  )
}

# createNumericInputWithPopover <- function(id, label, popoverContent, ...) {
#   tagList(
#     numericInput(id, label, ...),
#     bsPopover(id, title = label, content = popoverContent, trigger = "hover")
#   )
# }