plotGenotypequalityPerGenotype <- function(sce){
  genotype.matrix.filtered <- as.data.frame(t(assay(altExp(sce), 'Genotype')))
  genoqual.matrix.filtered <- as.data.frame(t(assay(altExp(sce), 'Genoqual')))
  colnames(genotype.matrix.filtered) <- paste0(rowData(altExp(sce))$Gene, ':', rowData(altExp(sce))$id)
  colnames(genoqual.matrix.filtered) <- paste0(rowData(altExp(sce))$Gene, ':', rowData(altExp(sce))$id)
  tmp.1 <- genotype.matrix.filtered %>% as.data.frame() %>% t() %>% melt(varnames = c('Variant', 'Cell'), value.name = 'Genotype')
  tmp.2 <-genoqual.matrix.filtered %>% as.data.frame() %>% t() %>% melt(varnames = c('Variant', 'Cell'), value.name = 'Genotype Quality')
  gt.df <- merge(tmp.1, tmp.2) %>%   
    dplyr::mutate(Genotype = ifelse(Genotype == 0, 'WT', 
                                    ifelse(Genotype == '1', 'Hom', 
                                           ifelse(Genotype == '2', 'Het', 'Missing'))))
  gt.df$Genotype <- factor(gt.df$Genotype, levels = c('Hom', 'Het', 'WT', 'Missing'))      
  
  var_plot5 <- ggplot(gt.df) +
    geom_violin(aes(x = Genotype, y =`Genotype Quality`, fill = Genotype), color = NA) +
    scale_fill_manual(values = mycols) +
    theme_default() +
    labs(title = 'Genotype Quality per Genotype (GATK)', x = NULL) +
    theme(panel.grid = element_blank()) +
    geom_hline(yintercept = 30, linetype = 'dashed')
  var_plot5
}

#})