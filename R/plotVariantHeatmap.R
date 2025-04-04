plotVariantHeatmap <- function(sce) {
  vaf.matrix.filtered = t(assay(altExp(sce, 'variants'), 'VAF'))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  
  genotype.matrix.filtered = t(assay(altExp(sce, 'variants'), 'Genotype'))
  source('./R/utils.R')
  
  # Create a color ramp for VAF
  colors_vaf <- circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292"))
  
  # Chromosome annotation
  column_ha <- HeatmapAnnotation(
    chr = factor(str_extract(colnames(vaf.matrix.filtered), 'chr(\\d|X|Y)+'), levels = chromosomes),
    col = list(chr = chr_palette)
  )
  
  # Process genotype matrix
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
  anno_bar <- anno_barplot(proportions, bar_width = 1, height = unit(3, "cm"),
                           gp = gpar(fill = c("#D44292", "#F6A97A", "#414487FF", "grey")))
  
  # Update column annotation with GT (%) annotation
  column_ha <- HeatmapAnnotation(
    chr = factor(str_extract(colnames(vaf.matrix.filtered), 'chr(\\d|X|Y)+'), levels = chromosomes),
    col = list(chr = chr_palette),
    'GT (%)' = anno_bar
  )
  
  # Draw the heatmap
  heatmap_plt <- Heatmap(matrix = vaf.matrix.filtered, 
                         name = 'VAF', 
                         col = colors_vaf,
                         show_column_dend = TRUE,
                         show_row_dend = FALSE, 
                         column_title = 'Filtered Variants',
                         row_title = 'Cells',  
                         top_annotation = column_ha)
  
  draw(heatmap_plt)
}
