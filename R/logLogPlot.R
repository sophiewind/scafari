#' Plot a log-log plot for read counts
#' 
#' This function generates a log-log plot of read counts using the data from a `SingleCellExperiment` object.
#'
#' @param sce A SingleCellExperiment object that contains read count data to be plotted. The read counts are extracted from an associated h5 file or similar data structure within the object.
#'
#' @return A ggplot object representing the log-log plot of read counts.

logLogPlot <- function(sce) {
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check if the counts assay is available
  if (!"counts" %in% assayNames(sce)) {
    stop("The SingleCellExperiment object must contain a 'counts' assay.")
  }
  
  read.counts <- as.data.frame(sce@assays@data$counts)
  # Ensure that the input is valid
  if (!is.data.frame(read.counts) || nrow(read.counts) == 0 || ncol(read.counts) == 0) {
    stop("The input must be a data frame.")
  }
  
  counts.per.cell <- read.counts %>% colSums()
  
  seq_plot3 <- counts.per.cell %>% as.data.frame() %>%
    rownames_to_column('barcode') %>%
    ggplot(data = .) +
    geom_point(aes(x = reorder(barcode, .),y  = .),stat = 'identity') +
    labs(x = 'barcodes', y = 'Number of reads', title = 'Reads per Barcode') +
    theme_default() +
    theme(axis.line = element_line(),
          axis.ticks.y = element_line(),
          panel.grid = element_blank(),
          axis.text.x = element_blank()) +
    scale_y_log10() +
    scale_x_discrete(limits=rev)
  return(seq_plot3)
}
