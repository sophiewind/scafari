plotGenotypeDistributionPie <- function(sce){
  genotype.matrix.filtered <- assay(altExp(sce, 'variants'), 'Genotype')
  genotype.matrix.filtered %>%
    table() %>%
    as.data.frame() %>%
    rename(Genotype = '.') %>%  # Renaming the first column to 'Genotype'
    dplyr::mutate(
      Genotype = factor(dplyr::case_when(
        Genotype == 0 ~ 'WT',
        Genotype == '1' ~ 'Hom',
        Genotype == '2' ~ 'Het',
        TRUE ~ 'Missing'),
        levels = c('Hom', 'Het', 'WT', 'Missing'))) %>%
    mutate(prop = round((Freq / sum(Freq) * 100))) %>%
    dplyr::arrange(desc(Genotype)) %>%  # Sort Genotype in descending order for correct pie order
    mutate(cumulative = cumsum(prop), 
           lab.ypos = cumulative - 0.5 * prop) %>% 
    ggplot(aes(x = 2, y = prop, fill = Genotype)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y", start = 0) +
    geom_text(aes(y = lab.ypos, label = paste0(prop, '%')), 
              color = "white", size = 5) +
    xlim(0.5, 2.5)   +
    scale_fill_manual(values = mycols) +
    theme_default() +
    xlim(0.5, 2.5) +
    labs(title = 'Genotype Distribution (%)') +
    theme(axis.text = element_blank(),
          axis.title = element_blank())
}
