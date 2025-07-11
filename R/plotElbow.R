#' This function takes a SingleCellExperiment object and variants of interest as
#'  input and plots an elbow plot to perform k-means later.
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#' @param variants.of.interest A vector specifying the variants of interest.
#'
#' @return ggplot object with elbow plot.
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp().
##' sce_filtered <- readRDS(system.file("extdata", "sce_filtered_demo.rds",
##' package = "scafari"))
#' plotElbow(
#'     sce = sce_filtered,
#'     variants.of.interest = c(
#'         "FLT3:chr13:28610183:A/G",
#'         "KIT:chr4:55599436:T/C",
#'         "TP53:chr17:7577427:G/A",
#'         "TET2:chr4:106158216:G/A"
#'     )
#' )
#' @export
plotElbow <- function(sce, variants.of.interest) {
    # Check if the input is a SingleCellExperiment object
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("The input must be a SingleCellExperiment object.")
    }

    # Check if variants.of.interest is a vector
    if (!is.vector(variants.of.interest)) {
        stop("The 'variants.of.interest' parameter must be a vector.")
    }

    # Extract the alternate experiment once
    alt_exp_data <- altExp(sce, "variants")

    # Attempt to extract and name the VAF matrix data
    vaf_data <- as.data.frame(t(assay(alt_exp_data, "VAF")))
    colnames(vaf_data) <- paste0(
        rowData(alt_exp_data)$Gene, ":",
        rowData(alt_exp_data)$id
    )
    vaf.matrix.filtered <- vaf_data

    # Check if the specified variants are present in the column names
    missing_variants <- setdiff(
        variants.of.interest,
        colnames(vaf.matrix.filtered)
    )
    if (length(missing_variants) > 0) {
        stop(
            "The following variants of interest are not present in the data: ",
            paste(missing_variants, collapse = ", ")
        )
    }

    # Filter, clean, and standardize the VAF matrix
    df <- na.omit(vaf.matrix.filtered[, variants.of.interest])
    df <- scale(df)

    # Determine the optimal number of clusters using the Elbow method
    plot <- fviz_nbclust(df, kmeans, method = "wss")

    return(plot)
}
