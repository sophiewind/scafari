#' Function: clusterVariantSelection
#' -------------------------------
#' This function takes selected variants and performs clustering on them.
#'
#' @param sce A SingleCellExperiment object containing the single-cell data on which clustering will be performed.
#' @param variants.of.interest A vector or list specifying the variants of interest to be selected for clustering.
#' @param n.clust An integer specifying the number of clusters.
#'
#' @return A list with k-means results and a ggplot-object.
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp()
#' clusterplot <- clusterVariantSelection(
#'   sce = sce_filtered,
#'   variants.of.interest = c(
#'     "FLT3:chr13:28610183:A/G",
#'     "KIT:chr4:55599436:T/C",
#'     "TP53:chr17:7577427:G/A",
#'     "TET2:chr4:106158216:G/A"
#'   ),
#'   n.clust = 4
#' )
#'
#' @export
clusterVariantSelection <- function(sce, variants.of.interest, n.clust) {
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }

  # Check that 'variants' altExp exists and contains a 'VAF' assay
  if (!"variants" %in% altExpNames(sce)) {
    stop("The SingleCellExperiment object must contain 
         'variants' as an alternate experiment.")
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

  df <- vaf.matrix.filtered[, variants.of.interest] # selected_variants()]
  df <- na.omit(df)
  df <- scale(df)

  # Determining Optimal Clusters
  kmeans_result <- kmeans(df, centers = n.clust, nstart = 25)


  # Generate the cluster plot and store it in shared_data
  gg.clust <- fviz_cluster(kmeans_result, data = df, ellipse.type = "norm", 
                           geom = "point") +
    theme_default()

  # Print the cluster plot for debugging
  print(gg.clust) # Display the cluster plot

  return(list(k_means = kmeans_result, clusterplot = gg.clust))
}
