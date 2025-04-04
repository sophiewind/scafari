plotClusterVAF <- function(sce, variants.of.interest, gg.clust){
  vaf.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  vaf.matrix.filtered <- vaf.matrix.filtered[,variants.of.interest] 
  
  #colnames(vaf.matrix.filtered) <- variants.of.interest 
  
  # add cluster information
  vaf.matrix.filtered.tmp <- vaf.matrix.filtered
  
  colnames(vaf.matrix.filtered.tmp) <-variants.of.interest 
  vaf.matrix.filtered.tmp$cluster <- paste0('c', gg.clust$data$cluster)
  vaf.matrix.filtered.tmp <- vaf.matrix.filtered.tmp %>%
    tidyr::pivot_longer(
      cols = c(-cluster),
      names_to = "variable",
      values_to = "value"
    ) %>% as.data.frame()
  vaf.matrix.filtered.tmp$variable <- factor(vaf.matrix.filtered.tmp$variable, levels = sort(variants.of.interest))
  
  p <- vaf.matrix.filtered.tmp %>%
    ggplot() +
    geom_violin(aes(x = cluster, y = value, fill = cluster), alpha = 0.5, col = NA) +
    geom_jitter(aes(x = cluster, y = value, col = cluster), size = 1) +
    labs(y = 'VAF', x = 'cluster')+
    facet_grid(~variable) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      title = element_text(size = 20),
      text = element_text(size = 18)
    )
  p
  return(p)
}
