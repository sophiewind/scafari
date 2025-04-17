#' Plot normalized read counts
#'
#' This function plots the normalize read counts in a bar chart. 
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#' 
#' @return ggplot object visualizing the normalized read counts in a bar chart.
#' 
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with 'counts' assay. 
#' plotNormalizedReadCounts(sce)
#'}
#'
#' @export 
plotNormalizedReadCounts <- function(sce) {
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check if 'normalized.counts' assay is available
  if (!"normalized.counts" %in% assayNames(sce)) {
    stop("The SingleCellExperiment object must contain a 'normalized.counts' assay.")
  }
  
  # Extract normalized read counts and verify it is not empty
  normalized.read.counts <- sce@assays@data$normalized.counts
  if (is.null(normalized.read.counts) || length(normalized.read.counts) == 0) {
    stop("The 'normalized.counts' assay is empty, cannot plot distribution.")
  }
  
  # Process the read counts and handle potential transformation issues
  tryCatch({
    # Compute normalized mean read counts per amplicon
    normalized.read.counts <- as.data.frame(rowMeans(normalized.read.counts)) %>%
      rownames_to_column('Amplicon')
    colnames(normalized.read.counts)[2] <- 'Normalized mean read\ncounts per amplicon'
    
    # Assign amplicon IDs from rowData
    if (!"id" %in% names(rowData(sce))) {
      stop("The rowData of SingleCellExperiment must contain an 'id' column for amplicon identification.")
    }
    normalized.read.counts$Amplicon <- rowData(sce)$id
    
    # Estimate font size for the plot
    font_size <- max(3, min(18, 3 + (ncol(normalized.read.counts) - 1) * (2 / (3000 - 1))))
    
    # Create the plot
    plot <- normalized.read.counts %>%
      ggplot() +
      geom_bar(aes(x = reorder(Amplicon, -`Normalized mean read\ncounts per amplicon`), y = `Normalized mean read\ncounts per amplicon`), stat = 'identity') +
      labs(x = '') +
      theme_default() +
      theme(axis.text.x = element_text(angle = -270, size = 10),
            axis.line = element_line())
    
    return(plot)
  }, error = function(e) {
    stop("An error occurred while transforming or plotting the data: ", e$message)
  })
}
