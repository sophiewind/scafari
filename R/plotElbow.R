plotElbow <- function(sce, variants.of.interest){
  # TODO schauen, ob ids anders als voi
  vaf.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  df <- vaf.matrix.filtered[,variants.of.interest] 
  df <- na.omit(df)
  df <- scale(df)
  
  # Determining Optimal Clusters - Elbow
  set.seed(123)
  plot <- fviz_nbclust(df, kmeans, method = "wss")
  
  # Store the plot in the reactive variable
  return(plot)
}
