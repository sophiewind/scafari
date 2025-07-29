#'
#' This function normalizes the read counts contained within a
#' `SingleCellExperiment` object.
#'
#' @param sce A SingleCellExperiment object that includes the assay data with
#' read counts to be normalized. The metadata within the object may also be
#' utilized for normalization purposes.
#'
#' @return SingleCellExperiment object  with normalized read counts.
#'
#' @references https://missionbio.github.io/mosaic/,
#' https://github.com/rachelgriffard/optima
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with 'counts' assay.
#' h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
#' h5 <- h5ToSce(h5_file_path)
#' sce <- h5$sce_amp
#' sce <- normalizeReadCounts(sce)
#'
#' @export
normalizeReadCounts <- function(sce) {
    checkSce(sce)

    # Check if the 'counts' assay is available and accessible
    if (!"counts" %in% assayNames(sce)) {
        stop("The SingleCellExperiment object must contain a 'counts' assay.")
    }

    # Extract read counts and metadata
    read.counts <- t(sce@assays@data$counts)
    metadata <- sce@metadata

    # Verify that metadata contains the expected element 'n_amplicons'
    if (is.null(metadata[["n_amplicons"]])) {
        stop("The metadata must contain 'n_amplicons'.")
    }

    # Check for proper dimension agreement
    if (nrow(read.counts) == 0 || ncol(read.counts) == 0) {
        stop("The 'counts' data must be non-empty.")
    }

    # Calculate row sum threshold
    if (nrow(read.counts) < 10) {
        stop("The 'counts' data must have at least 10 entries to compute a valid
        threshold.")
    }

    # Calculate row sum threshold
    rowsum_threshold <- sort(rowSums(read.counts), decreasing = TRUE)[10] / 10

    # Extract cells that will not be removed due to small total read counts
    keep_mask <- rowSums(read.counts) > rowsum_threshold |
        rowSums(read.counts) > 40 * as.numeric(metadata[["n_amplicons"]])

    # Normalize read counts
    read_counts_df_tmp <- read.counts / (apply(read.counts, 1, mean) + 1)

    # Normalize over cells using median reads per amp in cells kept
    read_counts_df_tmp <- t(t(read_counts_df_tmp) /
        (apply(
            read_counts_df_tmp[keep_mask, ], 2,
            median
        ) + 0.05))

    # Scale the normalized counts by a factor of 2
    read_counts_df_tmp <- read_counts_df_tmp * 2


    assays(sce, withDimnames = FALSE)$normalized.counts <- t(read_counts_df_tmp)

    # Convert back to data frame and return
    return(sce)
}
