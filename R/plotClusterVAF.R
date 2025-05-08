#' plot Cluster VAF
#' This function generates a plot to visualize variant allele frequency (VAF) in clusters based on selected variants of interest.
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#' @param variants.of.interest A vector specifying the variants of interest.
#' @param gg.clust An object containing clustering information.
#'
#' @return A ggplot object that visually represents the VAF in the clusters.
#' 
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp() and clusterplot is the output of clusterVariantSleection().
#' sce_filtered <- readRDS(system.file("extdata", "sce_filtered.rds", package = "scafari"))
#' clusterplot <- readRDS(system.file("extdata", "clusterplot.rds", package = "scafari"))
#' plotClusterVAF(
#'   sce = sce_filtered,
#'   variants.of.interest = c(
#'     "FLT3:chr13:28610183:A/G",
#'     "KIT:chr4:55599436:T/C",
#'     "TP53:chr17:7577427:G/A",
#'     "TET2:chr4:106158216:G/A"
#'   ),
#'   gg.clust = clusterplot$clusterplot
#' )
#'
#' @export 
plotClusterVAF <- function(sce, variants.of.interest, gg.clust){
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check that 'variants' altExp exists and contains a 'VAF' assay
  if (!"variants" %in% altExpNames(sce)) {
    stop("The SingleCellExperiment object must contain 'variants' as an alternate experiment.")
  }
  if (!"VAF" %in% assayNames(altExp(sce, "variants"))) {
    stop("The 'variants' alternate experiment must contain a 'VAF' assay.")
  }
  
  # Check that variants.of.interest is non-empty and exists in the data
  if (length(variants.of.interest) == 0) {
    stop("variants.of.interest must be a non-empty vector.")
  }
  
  vaf.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  
  if (!all(variants.of.interest %in% colnames(vaf.matrix.filtered))) {
    stop("All variants.of.interest must exist in the VAF matrix columns.")
  }
  
  
  vaf.matrix.filtered <- vaf.matrix.filtered[,variants.of.interest] 
  # Check if gg.clust contains necessary data and if it has correct dimensions
  if (!is.data.frame(gg.clust$data) || !"cluster" %in% names(gg.clust$data)) {
    stop("gg.clust$data must be a non-empty data frame with a 'cluster' column.")
  }
  if (nrow(gg.clust$data) != nrow(vaf.matrix.filtered)) {
    stop("The number of rows in gg.clust$data must match the number of rows in the VAF matrix.")
  }
  
  # add cluster information
  vaf.matrix.filtered.tmp <- vaf.matrix.filtered
  
  colnames(vaf.matrix.filtered.tmp) <-variants.of.interest 
  vaf.matrix.filtered.tmp$cluster <- paste0('c', gg.clust$data$cluster)
  vaf.matrix.filtered.tmp <- vaf.matrix.filtered.tmp %>%
    tidyr::pivot_longer(
      cols = c(-cluster),
      names_to = "variable",
      values_to = "value"
    ) %>% as.data.frame()
  vaf.matrix.filtered.tmp$variable <- factor(vaf.matrix.filtered.tmp$variable, levels = sort(variants.of.interest))
  
  p <- vaf.matrix.filtered.tmp %>%
    ggplot() +
    geom_violin(aes(x = cluster, y = value, fill = cluster), alpha = 0.5, col = NA) +
    geom_jitter(aes(x = cluster, y = value, col = cluster), size = 1) +
    labs(y = 'VAF', x = 'cluster')+
    facet_grid(~variable) +
    theme_default()
  p
  return(p)
}
