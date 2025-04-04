#' Plot panel uniformity
#' 
#' @param read.counts.norm dataframe with normalized read counts
#' @param amplicons amplicon information from h5 file
#'
#' @return Panel uniformity plot as ggplot object
plotPanelUniformity <- function(sce) {
  read.counts.norm <- as.data.frame(t(sce@assays@data$normalized.counts))
  amplicons <- as.data.frame(rowData(sce))
  # Ensure the input data is valid
  if (!is.data.frame(read.counts.norm) || !is.data.frame(amplicons)) {
    stop("Both inputs must be data frames.")
  }
  if (!all(c("id") %in% names(amplicons))) {
    stop("gene_anno_df must have an 'id' column.")
  }
  
  # Set column names based on gene annotation IDs
  colnames(read.counts.norm) <- amplicons$id
  
  # Melt the data frame for ggplot2
  read.counts.norm.melt <- melt(read.counts.norm)
  
  # Order the columns based on the mean of read counts
  tmp_m <- colMeans(read.counts.norm)
  tmp_m_o <- order(tmp_m)
  read.counts.norm.melt$variable <- factor(
    read.counts.norm.melt$variable, 
    levels = names(tmp_m[tmp_m_o])
  )
  
  # Calculate mean coverage
  mean_tmp <- mean(as.matrix(read.counts.norm))
  
  # Identify low uniformity amplicons
  low_uniformity_amps <- names(tmp_m[tmp_m < 0.2 * mean_tmp])
  
  # Add uniformity info for coloring
  read.counts.norm.melt <- read.counts.norm.melt %>%
    mutate(coverage = ifelse(variable %in% low_uniformity_amps, FALSE, TRUE))
  
  # Generate the ggplot
  plot <- ggplot(read.counts.norm.melt, aes(x = value, y = variable, color = coverage)) +
    geom_point() +
    scale_color_manual(values = c('TRUE' = '#414487FF', 'FALSE' = '#F66D7A')) +
    labs(x = 'Normalized read counts', y = 'Amplicons') +
    scale_x_continuous(limits = c(0, max(read.counts.norm.melt$value)), expand = c(0,0)) +
    theme_default() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line = element_line(),
          legend.position = 'none')
  
  return(plot)
}
