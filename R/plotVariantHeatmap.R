#' Plot Variant Heatmap
#' 
#' This function generates a heatmap to visualize variant data within a 
#' `SingleCellExperiment` object. The heatmap provides insights into the 
#' distribution and frequency of variants across different samples or cells.
#'
#' @param sce A `SingleCellExperiment` object containing single-cell variant data. 
#' The object should include an assay that holds variant information suitable 
#' for visualization in a heatmap.
#'
#' @return A `ggplot` object representing the heatmap of variant frequencies or 
#' presence across cells or groups.
#'
#' @examples
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with an appropriate variant assay.
#' variant_heatmap <- plotVariantHeatmap(sce)
#' print(variant_heatmap)
#' }
#'
#' @export 
plotVariantHeatmap <- function(sce) {
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check for the presence of 'variants' altExp and required assays
  if (!"variants" %in% altExpNames(sce)) {
    stop("The SingleCellExperiment object must contain 'variants' as an alternate experiment.")
  }
  if (!all(c("VAF", "Genotype") %in% assayNames(altExp(sce, "variants")))) {
    stop("The 'variants' alternate experiment must contain 'VAF' and 'Genotype' assays.")
  }
  
  # Extract and verify VAF and genotype matrices
  vaf.matrix.filtered <- t(assay(altExp(sce, 'variants'), 'VAF'))
  genotype.matrix.filtered <- t(assay(altExp(sce, 'variants'), 'Genotype'))
  
  if (nrow(vaf.matrix.filtered) == 0 || ncol(vaf.matrix.filtered) == 0) {
    stop("The VAF matrix is empty, cannot plot heatmap.")
  }
  if (nrow(genotype.matrix.filtered) == 0 || ncol(genotype.matrix.filtered) == 0) {
    stop("The Genotype matrix is empty, cannot plot heatmap.")
  }
  
  if (!'annotated' %in% names(metadata(altExp(sce_filtered)))){
    stop('The variants are not annotated. Please do this using `annotateVariants()`')
  }
  
  # Verify the presence of gene and id in rowData
  row_data <- rowData(altExp(sce, 'variants'))
  if (!all(c("Gene", "id") %in% names(row_data))) {
    stop("The rowData of 'variants' must contain 'Gene' and 'id' columns.")
  }
  
  colnames(vaf.matrix.filtered) <- paste0(row_data$Gene, ':', as.character(row_data$id))

  
  # Chromosome annotation
  column_ha <- HeatmapAnnotation(
    chr = factor(stringr::str_extract(colnames(vaf.matrix.filtered), 'chr(\\d|X|Y)+'), levels = chromosomes),
    col = list(chr = chr_palette)
  )
  
  # Process genotype matrix and check for expected columns
  gt_anno <- data.frame(WT = integer(), Het = integer(), Hom = integer(), Missing = integer())
  for (col in 1:ncol(genotype.matrix.filtered)) {
    wt <- sum(genotype.matrix.filtered[, col] == 0)
    het <- sum(genotype.matrix.filtered[, col] == 1)
    hom <- sum(genotype.matrix.filtered[, col] == 2)
    mis <- sum(genotype.matrix.filtered[, col] == 3)
    gt_anno[col, ] <- c(hom, het, wt, mis)
  }
  gt_anno$Total <- rowSums(gt_anno)
  
  proportions <- gt_anno %>%
    mutate(across(c(Hom, Het, WT, Missing), ~ . / Total * 100)) %>%
    dplyr::select(-Total)
  rownames(proportions) <- colnames(vaf.matrix.filtered)
  
  # Create a bar plot annotation for the genotype proportions
  anno_bar <- ComplexHeatmap::anno_barplot(proportions, bar_width = 1, height = unit(3, "cm"),
                                           gp = gpar(fill = c("#D44292", "#F6A97A", "#414487FF", "grey")))
  
  # Update column annotation with GT (%) annotation
  column_ha <- HeatmapAnnotation(
    chr = factor(stringr::str_extract(colnames(vaf.matrix.filtered), 'chr(\\d|X|Y)+'), levels = chromosomes),
    col = list(chr = chr_palette),
    'GT (%)' = anno_bar
  )
  
  # Draw the heatmap
  heatmap_plt <- ComplexHeatmap::Heatmap(matrix = vaf.matrix.filtered, 
                                         name = 'VAF', 
                                         col = colors_vaf,
                                         show_column_dend = TRUE,
                                         show_row_dend = FALSE, 
                                         column_title = 'Filtered Variants',
                                         row_title = 'Cells',  
                                         top_annotation = column_ha)
  
  ComplexHeatmap::draw(heatmap_plt)
}
