#' Plot Genotype Clusters
#' 
#' This function generates a plot to visualize genotype in clusters based on selected variants of interest.
#'
#' @param sce A SingleCellExperiment object containing the relevant data.
#' @param variants.of.interest A vector specifying the variants of interest.
#' @param gg.clust An object containing clustering information.
#'
#' @return A ggplot object that visually represents the clustering of genotypes based on the specified variants and clustering information.
#' 
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp() and clusterplot is the output of clusterVariantSleection().
#' plotClusterGenotype(sce = sce_filtered, 
#'                variants.of.interest = c("FLT3:chr13:28610183:A/G",
#'                                         "KIT:chr4:55599436:T/C",
#'                                         "TP53:chr17:7577427:G/A",
#'                                         "TET2:chr4:106158216:G/A"), 
#'                gg.clust = clusterplot$clusterplot)
#'
#' @export 
plotClusterGenotype <- function(sce, variants.of.interest, gg.clust){
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check that 'variants' altExp exists and contains a 'VAF' assay
  if (!"Genotype" %in% assayNames(altExp(sce, "variants"))) {
    stop("The 'variants' alternate experiment must contain a 'VAF' assay.")
  }
  
  # Check that variants.of.interest is non-empty and exists in the data
  if (length(variants.of.interest) == 0) {
    stop("variants.of.interest must be a non-empty vector.")
  }
  
  genotype.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'Genotype')))
  colnames(genotype.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  
  if (!all(variants.of.interest %in% colnames(genotype.matrix.filtered))) {
    stop("All variants.of.interest must exist in the VAF matrix columns.")
  }
  
  genotype.matrix.filtered <- genotype.matrix.filtered[,variants.of.interest] 
  
  
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
            theme_default()
  return(p)
}
