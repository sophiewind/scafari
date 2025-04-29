#' Plot genotype distribution pie
#' 
#' This function generates a plot to visualize the genotype distribution.
#' 
#' @param sce A SingleCellExperiment object containing the relevant data.
#' 
#' @return ggplot object visualizing the genotype distribution in a pie chart.
#' 
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp()
#' pie_chart <- plotGenotypeDistributionPie(sce)
#' print(pie_chart)
#'}
#'
#' @export
plotGenotypeDistributionPie <- function(sce) {
 # browser()
  mycols <- c('WT' = "#414487FF", 'Hom' = "#D44292", 'Het' = "#F6A97A",
              'Missing' = "#868686FF")
  
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check if the 'variants' altExp exists and contains a 'Genotype' assay
  if (!"variants" %in% altExpNames(sce)) {
    stop("The SingleCellExperiment object must contain 'variants' as an alternate experiment.")
  }
  if (!"Genotype" %in% assayNames(altExp(sce, "variants"))) {
    stop("The 'variants' alternate experiment must contain a 'Genotype' assay.")
  }
  
  # Extract the genotype matrix and check if it is empty
  genotype.matrix.filtered <- as.data.frame(assay(altExp(sce, 'variants'), 'Genotype'))
  if (is.null(genotype.matrix.filtered) || length(genotype.matrix.filtered) == 0) {
    stop("The Genotype matrix is empty, cannot plot distribution.")
  }
  
  
  # Transform and plot the genotype distribution as a pie chart
  #tryCatch({
    test <- genotype.matrix.filtered %>%
            table() %>%
            rename(Genotype = '.')
    test.2 <- test %>%  # Renaming the first column to 'Genotype'
            dplyr::mutate(
              Genotype = factor(dplyr::case_when(
                Genotype == 0 ~ 'WT',
                Genotype == '1' ~ 'Hom',
                Genotype == '2' ~ 'Het',
                TRUE ~ 'Missing'),
                levels = c('Hom', 'Het', 'WT', 'Missing'))) 
    test3 <- test.2 %>%
            mutate(prop = round((Freq / sum(Freq) * 100))) %>%
            dplyr::arrange(desc(Genotype)) %>%
            mutate(cumulative = cumsum(prop), 
                   lab.ypos = cumulative - 0.5 * prop)
    test3
    # genotype.matrix.filtered %>%
    #   table() %>%
    #   as.data.frame() %>%
    #   rename(Genotype = '.') %>%  # Renaming the first column to 'Genotype'
    #   dplyr::mutate(
    #     Genotype = factor(dplyr::case_when(
    #       Genotype == 0 ~ 'WT',
    #       Genotype == '1' ~ 'Hom',
    #       Genotype == '2' ~ 'Het',
    #       TRUE ~ 'Missing'),
    #       levels = c('Hom', 'Het', 'WT', 'Missing'))) %>%
    #   mutate(prop = round((Freq / sum(Freq) * 100))) %>%
    #   dplyr::arrange(desc(Genotype)) %>%
    #   mutate(cumulative = cumsum(prop), 
    #          lab.ypos = cumulative - 0.5 * prop) %>% 
    #   ggplot(aes(x = 2, y = prop, fill = Genotype)) +
    #   geom_bar(stat = "identity", width = 1) +
    #   coord_polar(theta = "y", start = 0) +
    #   geom_text(aes(y = lab.ypos, label = paste0(prop, '%')), 
    #             color = "white", size = 5) +
    #   xlim(0.5, 2.5) +
    #   scale_fill_manual(values = mycols) +
    #   theme_default() +
    #   xlim(0.5, 2.5) +
    #   labs(title = 'Genotype Distribution (%)') +
    #   theme(axis.text = element_blank(),
    #         axis.title = element_blank())
  # }, error = function(e) {
  #   stop("An error occurred while transforming or plotting the data: ", e$message)
  # })
}
