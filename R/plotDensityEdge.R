plotDensityEdge <- function(sce, variants.of.interest, min.pts) {
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
  
  vaf.matrix.filtered <- as.data.frame(t(assay(altExp(sce, "variants"), "VAF")))
  colnames(vaf.matrix.filtered) <- 
    paste0(rowData(altExp(sce, "variants"))$Gene, ":", 
           rowData(altExp(sce, "variants"))$id)
  
  if (!all(variants.of.interest %in% colnames(vaf.matrix.filtered))) {
    stop("All variants.of.interest must exist in the VAF matrix columns.")
  }
  
  df <- vaf.matrix.filtered[, variants.of.interest]
  df <- na.omit(df)
  df <- scale(df)
  
  dbscan::kNNdistplot(df, k = min.pts)
}
