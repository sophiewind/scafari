#' Plot Genotype Quality per Genotype
#' 
#' This function generates a plot to assess the quality of genotypes within a 
#' `SingleCellExperiment` object. It uses the 'genotype_quality' assay or any 
#' relevant assay to visualize the distribution of genotype quality metrics 
#' across different genotypes.
#'
#' @param sce A `SingleCellExperiment` object containing the single-cell data. 
#' This object should include a 'genotype_quality' assay or similar data which 
#' provides quality metrics for each genotype.
#'
#' @return A `ggplot` object visualizing the distribution of genotype quality 
#' across different genotypes.
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with 'genotype_quality' assay.
#' sce_filtered <- readRDS(system.file("extdata", "sce_filtered.rds", package = "scafari"))
#' genotype_quality_plot <- plotGenotypequalityPerGenotype(sce_filtered)
#' print(genotype_quality_plot)
#'
#' @export 
plotGenotypequalityPerGenotype <- function(sce) {
  tryCatch({
    # Check if the 'sce' object is valid
    if (is.null(sce) || !(inherits(sce, "SingleCellExperiment"))) {
      stop("Input 'sce' is invalid. It must be a SingleCellExperiment object.")
    }
    
    # Retrieve genotype and genotype quality matrices
    genotype.matrix.filtered <- tryCatch({
      as.data.frame(t(assay(altExp(sce), 'Genotype')))
    }, error = function(e) {
      stop("Failed to retrieve or transpose genotype matrix: ", e$message)
    })
    
    genoqual.matrix.filtered <- tryCatch({
      as.data.frame(t(assay(altExp(sce), 'Genoqual')))
    }, error = function(e) {
      stop("Failed to retrieve or transpose genotype quality matrix: ", e$message)
    })
    
    # Check if the matrices are non-empty
    if (nrow(genotype.matrix.filtered) == 0 || ncol(genotype.matrix.filtered) == 0 ||
        nrow(genoqual.matrix.filtered) == 0 || ncol(genoqual.matrix.filtered) == 0) {
      stop("One or both matrices (genotype/genotype quality) are empty.")
    }
    
    # Assign column names to matrices
    tryCatch({
      colnames(genotype.matrix.filtered) <- paste0(rowData(altExp(sce))$Gene, ':', rowData(altExp(sce))$id)
      colnames(genoqual.matrix.filtered) <- paste0(rowData(altExp(sce))$Gene, ':', rowData(altExp(sce))$id)
    }, error = function(e) {
      stop("Failed to assign column names to matrices: ", e$message)
    })
    
    # Data transformation steps
    tmp.1 <- tryCatch({
      genotype.matrix.filtered %>% as.data.frame() %>% t() %>% melt(varnames = c('Variant', 'Cell'), value.name = 'Genotype')
    }, error = function(e) {
      stop("Failed during melting genotype matrix: ", e$message)
    })
    
    tmp.2 <- tryCatch({
      genoqual.matrix.filtered %>% as.data.frame() %>% t() %>% melt(varnames = c('Variant', 'Cell'), value.name = 'Genotype Quality')
    }, error = function(e) {
      stop("Failed during melting genotype quality matrix: ", e$message)
    })
    
    # Merge and mutate data
    gt.df <- tryCatch({
      merge(tmp.1, tmp.2) %>%
        dplyr::mutate(Genotype = ifelse(Genotype == 0, 'WT', 
                                        ifelse(Genotype == '1', 'Hom', 
                                               ifelse(Genotype == '2', 'Het', 'Missing'))))
    }, error = function(e) {
      stop("Failed to merge and mutate data: ", e$message)
    })
    
    # Set Genotype as a factor
    gt.df$Genotype <- factor(gt.df$Genotype, levels = c('Hom', 'Het', 'WT', 'Missing'))
    
    # Plot creation
    var_plot5 <- tryCatch({
      ggplot(gt.df) +
        geom_violin(aes(x = Genotype, y = `Genotype Quality`, fill = Genotype), color = NA) +
        scale_fill_manual(values = mycols) +
        theme_default() +
        labs(title = 'Genotype Quality per Genotype (GATK)', x = NULL) +
        theme(panel.grid = element_blank()) +
        geom_hline(yintercept = 30, linetype = 'dashed')
    }, error = function(e) {
      stop("Failed to create plot: ", e$message)
    })
    
    return(var_plot5)
    
  }, error = function(e) {
    message("An error occurred: ", e$message)
    return(NULL)
  })
}