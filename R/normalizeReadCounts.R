#' Normalize read counts
#' 
#' @param read.counts datafrane with read counts from h5 file
#' @param metadata dataframe with metadata from h5 file
#'
#' @return normalized read counts
#' @references TODO scafari
normalizeReadCounts <- function(sce) {
  read.counts <- t(sce@assays@data$counts)
  metadata <- sce@metadata

  # Calculate row sum threshold
  rowsum_threshold <- sort(rowSums(read.counts), decreasing = TRUE)[10] / 10
  
  # Extract cells that will not be removed due to small total read counts
  keep_mask <- rowSums(read.counts) > rowsum_threshold | 
    rowSums(read.counts) > 40 * as.numeric(metadata[['n_amplicons']])
  
  # Normalize read counts
  read_counts_df_tmp <- read.counts / (apply(read.counts, 1, mean) + 1)
  
  # Normalize over cells using median reads per amp in cells kept
  read_counts_df_tmp <- t(t(read_counts_df_tmp) / (apply(read_counts_df_tmp[keep_mask, ], 2, median) + 0.05))
  
  # Scale the normalized counts by a factor of 2
  read_counts_df_tmp <- read_counts_df_tmp * 2
  
  # Convert back to data frame and return
  return(as.data.frame(t(read_counts_df_tmp)))
  
}
