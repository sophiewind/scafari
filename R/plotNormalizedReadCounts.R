#' Plot normalized read counts
#'
#' This function plots the normalize read counts in a bar chart. 
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#' 
#' @return ggplot object visualizing the normalized read counts in a bar chart.

plotNormalizedReadCounts <- function(sce) {
  normalized.read.counts <- sce@assays@data$normalized.counts
  
  # Normalized mean read counts per amplicon
  normalized.read.counts <- as.data.frame(normalized.read.counts %>% rowMeans()) %>% 
    rownames_to_column('Amplicon')
  colnames(normalized.read.counts)[2] <- 'Normalized mean read\ncounts per amplicon'
  
  normalized.read.counts$Amplicon <- rowData(sce)$id
  font_size <- max(3, min(18, 3 + (ncol(normalized.read.counts
  ) - 1) * (2 / (3000 - 1))))
  plot <- normalized.read.counts%>% ggplot(.) +
    geom_bar(aes(x = reorder(Amplicon, -`Normalized mean read\ncounts per amplicon`), y = `Normalized mean read\ncounts per amplicon`), stat = 'identity') +
    labs(x = '') +
    theme_default() +
    theme_default() +
    theme(axis.text.x = element_text(angle = -270, size = 10),
          axis.line = element_line())
  # plot
  return(plot)
}
