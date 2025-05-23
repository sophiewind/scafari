#' Plot Distribution of Amplicons
#' 
#' This function generates a plot to visualize the distribution of amplicons within a `SingleCellExperiment` object.
#'
#' @param sce A SingleCellExperiment object that contains the assay data, including amplicon information to be plotted.
#'
#' @return A ggplot object representing the distribution of amplicons.
#' 
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with 'counts' assay.
#' h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
#' h5 <- h5ToSce(h5_file_path)
#' sce <- h5$sce_amp 
#' plotAmpliconDistribution(sce)
#'
#' @export 
plotAmpliconDistribution <- function(sce) {
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  amps <- as.data.frame(rowData(sce))%>% 
    mutate(Gene = vapply(strsplit(id, '_'), function(x) x[3], character(1))) %>% 
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE) %>% 
    as.data.frame()
  
  # Ensure the input data is a data frame
  if (!is.data.frame(amps)) {
    stop("gene_anno_df must be a data frame.")
  }
  
  data(ideoCyto, package = "biovizBase")

  
  amps$tooltip <- paste("Gene:", amps$Gene, "<br>ID:", amps$id)
  
  plot <- ggbio::autoplot((ideoCyto$hg19), layout = "karyogram", cytobands = TRUE) +
    geom_segment(data = amps, aes(x = start, xend = start + width,
                                  y = -2, 
                                  yend = 12, 
                                  text = tooltip
    ), color = "blue", size = 1) +
    theme(panel.background = element_blank())
  
  return(plot)
}
