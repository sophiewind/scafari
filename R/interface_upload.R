createUploadUI <- function() {
    tabPanel(
        "Upload",
        fluidPage(
            width = 3,
            h3("Data Upload"),
            fileInput("upload", "Upload HDF5 File", accept = ".h5"),
            verbatimTextOutput("class_output"),
            verbatimTextOutput("file_contents"),
            div(id = "text"),
            tableOutput("files")
        )
    )
}
