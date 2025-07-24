#' Plot Panel Uniformity
#'
#' This function generates a plot to assess the uniformity of panel reads
#' in a `SingleCellExperiment` object. It uses read counts stored in the
#' 'counts' assay to visualize the distribution and variability of reads.
#'
#' @param sce A `SingleCellExperiment` object that contains single-cell data
#' with a 'counts' assay. This object should be pre-processed to ensure
#' the data is appropriate for uniformity analysis.
#' @param interactive in interactive mode an plotly is returned.
#' @return A `ggplot` object visualizing the uniformity of panel reads.
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with 'counts' assay.
#' h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
#' h5 <- h5ToSce(h5_file_path)
#' sce <- h5$sce_amp
#' sce <- normalizeReadCounts(sce = sce)
#' uniformity_plot <- plotPanelUniformity(sce)
#' print(uniformity_plot)
#'
#' @export
plotPanelUniformity <- function(sce, interactive = FALSE) {
    checkSce(sce)
    
    if (!inherits(interactive, "logical")) {
        stop("`interactive` must be logical.")}
    
    if (!"normalized.counts" %in% assayNames(sce)) {
        stop("The SingleCellExperiment object must contain a ",
                    "'normalized.counts' assay.")}

    # Extract and verify normalized read counts and amplicon data
    read.counts.norm <- as.data.frame(t(sce@assays@data$normalized.counts))
    amplicons <- as.data.frame(rowData(sce))

    # Ensure the assays and rowData contain the necessary data
    if (nrow(read.counts.norm) == 0 || ncol(read.counts.norm) == 0) {
        stop("The 'normalized.counts' assay is empty, cannot plot panel",
            "uniformity.")}
    if (nrow(amplicons) == 0 || !"id" %in% names(amplicons)) {
        stop("The rowData must contain non-empty data with an 'id' column.")}

    colnames(read.counts.norm) <- amplicons$id
    read.counts.norm.melt <- reshape2::melt(read.counts.norm)

    # Order the columns based on the mean of read counts
    tmp_m <- colMeans(read.counts.norm, na.rm = TRUE)
    tmp_m_o <- order(tmp_m)
    read.counts.norm.melt$variable <- factor(read.counts.norm.melt$variable,
        levels = names(tmp_m[tmp_m_o]))

    mean_tmp <- mean(as.matrix(read.counts.norm), na.rm = TRUE)

    # Identify low uniformity amplicons
    low_uniformity_amps <- names(tmp_m[tmp_m < 0.2 * mean_tmp])

    # Add uniformity info for coloring
    read.counts.norm.melt <- read.counts.norm.melt %>% mutate(coverage = ifelse(
        variable %in% low_uniformity_amps, FALSE, TRUE))

    plot <- ggplot(read.counts.norm.melt, aes(
        x = value, y = variable, color = coverage)) +
        geom_point() +
        scale_color_manual(values = c("TRUE" = "#414487FF", "FALSE" = "#F66D7A"
        )) +
        labs(x = "Normalized read counts", y = "Amplicons") +
        scale_x_continuous(limits = c(0, max(read.counts.norm.melt$value,
            na.rm = TRUE)), expand = c(0, 0)) +
        theme_default() +
        theme(panel.grid = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank(), axis.line = element_line(),
            legend.position = "none")

    if (interactive) {
        plot <- ggplotly(plot)
        return(plot)
    } else {
        return(plot)
    }
}
