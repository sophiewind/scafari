#' plot Cluster VAF Map
#' This function generates a plot to visualize variant allele frequency (VAF) in clusters based on selected variants of interest with clusters in the background.
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#' @param variants.of.interest A vector specifying the variants of interest.
#' @param gg.clust An object containing clustering information.
#'
#' @return A ggplot object that visually represents the VAF in the clusters with clusters in the background.
#' 
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp() and clusterplot is the output of clusterVariantSleection().
#' plotClusterVAFMap(sce = sce_filtered, 
#'                variants.of.interest = c("FLT3:chr13:28610183:A/G",
#'                                         "KIT:chr4:55599436:T/C",
#'                                         "TP53:chr17:7577427:G/A",
#'                                         "TET2:chr4:106158216:G/A"), 
#'                gg.clust = clusterplot$clusterplot)
#'}
#'
#' @export 
plotClusterVAFMap <- function(sce, variants.of.interest, gg.clust){
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
  
  # Check if gg.clust contains necessary data
  if (!is.data.frame(gg.clust$data) || nrow(gg.clust$data) == 0) {
    stop("gg.clust$data must be a non-empty data frame.")
  }
  
  # Ensure gg.clust contains x and y labels
  if (!all(c("x", "y") %in% names(gg.clust$labels))) {
    stop("gg.clust must contain 'x' and 'y' labels for plotting.")
  }
  
  
  vaf.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  rownames(vaf.matrix.filtered) <- paste0('cell', rownames(vaf.matrix.filtered)) # TODO cell ids
  
  if (!all(variants.of.interest %in% colnames(vaf.matrix.filtered))) {
    stop("All variants.of.interest must exist in the VAF matrix columns.")
  }
 
  vaf.matrix.filtered <- vaf.matrix.filtered[,variants.of.interest] 
  rownames(gg.clust$data) <- paste0('cell', rownames(gg.clust$data))
  
  merged <- merge(gg.clust$data, vaf.matrix.filtered, by = 0)
  col_names <- colnames(merged)
  gene_cols <- col_names[grepl("[^:]:chr[[:digit:]]+", col_names, perl = T)]
  
  merged.var <- merged %>%
    tidyr::pivot_longer(
      cols = all_of(gene_cols),
      names_to = "variant", 
      values_to = "VAF" 
    )
  merged.var$variant <- factor(merged.var$variant, levels = sort(variants.of.interest))
  
  # Define the color function with gradient and -1 as grey
  color_func <- function(value) {
    ifelse(value == -1, "#BEBEBE",  # Grey for -1
           circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292"))(value))}
  
  # Create a sequence of breakpoints covering the data range you want
  breakpoints <- seq(-1, 100, length.out = 101)  # Adjust based on your data, now includes -1
  
  # Generate the corresponding colors
  colors <- color_func(breakpoints)
  
  # legend
  data <- data.frame(
    value = c(-1, seq(0, 100, by = 25)),
    label = c("Missing", "0", "25", "50", "75", "100")
  )
  
  
  p <- ggplot(merged.var, aes(x = x, y = y, color = VAF)) +
      geom_polygon(
        aes(x = x, y = y, group = cluster, fill = cluster), 
        alpha = 0.3,
        stat = "ellipse", 
        type = "norm", 
        level = 0.95, 
        segments = 51,
        na.rm = FALSE) +
      geom_point() +
      scale_color_gradientn(colors = colors, values = scales::rescale(breakpoints)) +
      geom_point() +
      facet_grid(~factor(variant, levels = sort(variants.of.interest))) +
      theme_default()
  return(p)
}
