plotClusterGenotype <- function(sce, variants.of.interest, gg.clust){
  genotype.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'Genotype')))
  colnames(genotype.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  
    #genotype.matrix.filtered <<- genotype.matrix.filtered
  #variants.of.interest <- variants.of.interest
  #gg.clust <<- gg.clust
  genotype.matrix.filtered <- genotype.matrix.filtered[,variants.of.interest] 
  #colnames(genotype.matrix.filtered) <- variants.of.interest   # TODO important, too
  
  # add cluster information
  genotype.matrix.filtered.tmp <- as.data.frame(genotype.matrix.filtered)
  
  genotype.matrix.filtered.tmp$cluster <- paste0('c', gg.clust$data$cluster)
  genotype.matrix.filtered.tmp <- melt(genotype.matrix.filtered.tmp)
  
  
  gt <- genotype.matrix.filtered.tmp %>%
    #melt() %>%
    mutate(variable = factor(variable, levels = sort(variants.of.interest))) %>% 
    mutate(value = as.factor(value)) %>%
    mutate(Genotype = factor(dplyr::case_when(
      value == 0 ~ "WT",
      value == 1 ~ "Het",
      value == 2 ~ "Hom",
      TRUE ~ "Missing"  # Fallback for unexpected values
    ), levels = c('Hom', 'Het', 'WT', 'Missing')))
  
  recode_function <- function(x) {
    ifelse(x == 0, "WT",
           ifelse(x == 1, "Het",
                  ifelse(x == 2, "Hom",
                         ifelse(x == 3, "Missing", x))))
  }
  
  # Plot barplot
  p <- gt %>%
            ggplot() +
            geom_bar(aes(x = cluster, fill = Genotype), col = NA) +
            #geom_jitter(aes(x = cluster, y = value, col = cluster), size = 1) +
            #geom_boxplot(aes(x = cluster, y = value), outliers = F) +
            scale_fill_manual(values = mycols.ngt) +
            #theme_default() +
            #labs(y = 'VAF', x = variant.of.interest)+
            facet_grid(~variable) +
            theme(
              panel.grid = element_blank(),
              panel.background = element_blank(),
              title = element_text(size = 20),
              text = element_text(size = 18))
  return(p)
}
