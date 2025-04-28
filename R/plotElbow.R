#' This function takes a SingleCellExperiment object and variants of interest as input and plots an elbow plot to perform k-means later.
#' 
#' @param sce A SingleCellExperiment object containing the relevant data.
#' @param variants.of.interest A vector specifying the variants of interest.
#' 
#' @return ggplot object with elbow plot.
#' 
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp().
#' plotElbow(sce = sce_filtered, 
#'                variants.of.interest = c("FLT3:chr13:28610183:A/G",
#'                                         "KIT:chr4:55599436:T/C",
#'                                         "TP53:chr17:7577427:G/A",
#'                                         "TET2:chr4:106158216:G/A"))
#'}
#'
#' @export 
plotElbow <- function(sce, variants.of.interest){
  # Check if the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check if variants.of.interest is a vector
  if (!is.vector(variants.of.interest)) {
    stop("The 'variants.of.interest' parameter must be a vector.")
  }
  
  # Extract the VAF matrix
  tryCatch({
    vaf.matrix.filtered <- as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  }, error = function(e) {
    stop("Error extracting VAF data: ", e$message)
  })
  
  # Handle potential errors during column name assignment
  tryCatch({
    colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  }, error = function(e) {
    stop("Error assigning column names: ", e$message)
  })
  
  # Check if the specified variants are present in the column names
  missing_variants <- setdiff(variants.of.interest, colnames(vaf.matrix.filtered))
  if (length(missing_variants) > 0) {
    stop("The following variants of interest are not present in the data: ", paste(missing_variants, collapse = ", "))
  }
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
