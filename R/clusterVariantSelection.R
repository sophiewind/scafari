clusterVariantSelection <- function(sce, variants.of.interest, n.clust){
  vaf.matrix.filtered <- as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  df <- vaf.matrix.filtered[, variants.of.interest] #selected_variants()] 
  df <- na.omit(df)
  df <- scale(df)
  
  # Determining Optimal Clusters
  kmeans_result <- kmeans(df, centers = n.clust, nstart = 25) 
  

  # Generate the cluster plot and store it in shared_data
  gg.clust <- fviz_cluster(kmeans_result, data = df, ellipse.type = "norm", geom = 'point') + 
             theme_default()
  
  # Print the cluster plot for debugging
  print(gg.clust)  # Display the cluster plot
  
  return(list(k_means = kmeans_result, clusterplot = gg.clust))
}
