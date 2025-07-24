#' Plot normalized read counts
#'
#' This function plots the normalize read counts in a bar chart.
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#'
#' @return ggplot object visualizing the normalized read counts in a bar chart.
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with 'counts' assay.
#' h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
#' h5 <- h5ToSce(h5_file_path)
#' sce <- h5$sce_amp
#' sce <- normalizeReadCounts(sce = sce)
#' plotNormalizedReadCounts(sce)
#'
#' @export
plotNormalizedReadCounts <- function(sce) {
    checkSce(sce)

    # Check if 'normalized.counts' assay is available
    if (!"normalized.counts" %in% assayNames(sce)) {
        stop("The SingleCellExperiment object must contain a ",
            "'normalized.counts' assay.")}

    normalized.read.counts <- sce@assays@data$normalized.counts
    if (is.null(normalized.read.counts) ||
        length(normalized.read.counts) == 0) {
        stop("The 'normalized.counts' assay is empty,
            cannot plot distribution.")}

    # Compute normalized mean read counts per amplicon
    normalized.read.counts <-
        as.data.frame(rowMeans(normalized.read.counts)) %>%
        rownames_to_column("Amplicon")
    colnames(normalized.read.counts)[2] <-
        "Normalized_mean_read_counts_per_amplicon"

    if (!"id" %in% names(rowData(sce))) {
        stop("The rowData of SingleCellExperiment must contain an 'id' ",
                    "column for amplicon identification.")}
    normalized.read.counts$Amplicon <- rowData(sce)$id

    # Estimate font size for the plot
    font_size <- max(3, min(18, 3 + (ncol(normalized.read.counts) - 1) *
        (2 / (3000 - 1))))

    plot <- normalized.read.counts %>%
        ggplot() + geom_bar(aes(x = reorder(Amplicon,
                    -`Normalized_mean_read_counts_per_amplicon`),
                y = `Normalized_mean_read_counts_per_amplicon`),
            stat = "identity") +
        labs(x = "", y = "Normalized mean read\ncounts per amplicon") +
        theme_default() +
        theme(axis.text.x = element_text(angle = -270, size = 8),
            axis.line = element_line())
    return(plot)
}
