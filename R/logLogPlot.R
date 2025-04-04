#' Plot a log-log plot for read counts
#' 
#' @param read.counts dataframe with read counts from h5 file
#'
#' @return log-log plot as ggplot object
logLogPlot <- function(sce) {
  read.counts <- as.data.frame(sce@assays@data$counts)
  # Ensure that the input is valid
  if (!is.data.frame(read.counts)) {
    stop("The input must be a data frame.")
  }
  
  counts.per.cell <- read.counts %>% colSums()
  
  seq_plot3 <- counts.per.cell %>% as.data.frame() %>%
    rownames_to_column('barcode') %>%
    ggplot(data = .) +
    geom_point(aes(x = reorder(barcode, .),y  = .),stat = 'identity') +
    labs(x = 'barcodes', y = 'Number of reads', title = 'Reads per Barcode') +
    theme_default() +
    theme(axis.line = element_line(),
          axis.ticks.y = element_line(),
          panel.grid = element_blank(),
          axis.text.x = element_blank()) +
    scale_y_log10() +
    scale_x_discrete(limits=rev)
  
  seq_plot3
  
  # # Calculate counts per cell
  # counts_per_cell <- rowSums(read_counts_df)
  # 
  # # Create a data frame for plotting
  # plot_data <- data.frame(
  #   barcode = rownames(read_counts_df),
  #   counts = counts_per_cell
  # )
  # 
  # # Generate the plot
  # seq_plot3 <- ggplot(data = plot_data) +
  #   geom_point(aes(x = reorder(barcode, counts), y = counts), stat = 'identity') +
  #   labs(x = 'Barcodes', y = 'Number of Reads', title = 'Reads per Barcode') +
  #   theme_minimal() + # Using theme_minimal for simplicity, replace with theme_default() if defined
  #   theme(
  #     axis.text.x = element_blank(),
  #     axis.ticks.x = element_blank(),
  #     panel.grid = element_blank()
  #   ) +
  #   scale_y_log10() +
  #   scale_x_discrete(limits = rev(rownames(read_counts_df)))
  # 
  return(seq_plot3)
}
