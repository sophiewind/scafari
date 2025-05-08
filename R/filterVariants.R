#' Function: filterVariants
#' -------------------------------
#' This function takes a SingleCellExperiment object as input and performs variant filtering
#' @param depth.threshold A numeric value specifying the minimum read depth required.
#' @param genotype.quality.threshold A numeric value specifying the minimum genotype quality score.
#' @param vaf.ref A numeric value specifying the variant allele frequency threshold for wild-type alleles.
#' @param vaf.het A numeric value specifying the variant allele frequency threshold for heterozygous variants.
#' @param vaf.hom A numeric value specifying the variant allele frequency threshold for homozygous variants.
#' @param min.cell A numeric value indicating the minimum number of cells with a genotype other than "missing".
#' @param min.mut.cell A numeric value indicating the minimum number of mutated (genotype either "homozygous" or "heterozygous") cells.
#' @param se.var The SummarizedExperiment object containing variant data which will be filtered.
#' @param sce SingleCellExperiment object containing the single-cell data TODO.
#' @param shiny A logical flag indicating whether the function is being run in a Shiny application
#' context. Default is FALSE. 
#' 
#' @return A list containing the following elements:
#' \describe{
#'   \item{vaf.matrix.filtered}{Variant allele frequencies after filtering.}
#'   \item{genotype.matrix.filtered}{Genotype information after filtering.}
#'   \item{read.counts.df.norm.filtered}{Normalized read counts for variants retained after filtering.}
#'   \item{variant.ids.filtered}{A vector of the variant IDs that were retained after filtering.}
#'   \item{genoqual.matrix.filtered}{Genotype qualities for variants retained after filtering.}
#'   \item{cells.keep}{A vector of cell identifiers for those cells retained after filtering.}
#' } 
#' 
#' @examples
#' \dontrun{
#' filteres <- filterVariants(depth.threshold = 10,
#' genotype.quality.threshold = 30,
#' vaf.ref = 5, 
#' vaf.het =  35, 
#' vaf.hom = 95, 
#' min.cell = 50,
#' min.mut.cell = 1,
#' se.var = se.var,
#' sce = sce,
#' shiny = FALSE)
#' se.f <- SummarizedExperiment(assays = list(VAF = t(filteres$vaf.matrix.filtered), 
#'                                            Genotype = t(filteres$genotype.matrix.filtered),
#'                                            Genoqual = t(filteres$genoqual.matrix.filtered)),
#'                              rowData = filteres$variant.ids.filtered,
#'                              colData = filteres$cells.keep)
#' 
#' # Filter out cells in sce object
#' # Find the indices of the columns to keep
#' indices_to_keep <- match(filteres$cells.keep, SummarizedExperiment::colData(sce)[[1]], nomatch = 0)
#' 
#' # Subset the SCE using these indices
#' sce_filtered <- sce[, indices_to_keep]
#' }
#' 
#' @references https://missionbio.github.io/mosaic/, https://github.com/rachelgriffard/optima
filterVariants <- function(depth.threshold = numeric(),
                           genotype.quality.threshold = numeric(),
                           vaf.ref = numeric(),
                           vaf.het = numeric(),
                           vaf.hom = numeric(),
                           min.cell = numeric(),
                           min.mut.cell = numeric(),
                           se.var,
                           sce,
                           shiny = FALSE){
  message('start variant filtering...\n')
  # Validate SingleCellExperiment inputs
  if (!inherits(se.var, "SummarizedExperiment")) {
    stop("se.var must be a SummarizedExperiment object.")
  }
  
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("sce must be a SingleCellExperiment object.")
  }
  
  # Check for required assays in se.var
  required_assays <- c("Depth", "VAF", "Genoqual", "Genotype")
  missing_assays <- required_assays[!(required_assays %in% assayNames(se.var))]
  if (length(missing_assays) > 0) {
    stop("Missing required assays in se.var: ", paste(missing_assays, collapse = ", "))
  }
  
  # Check for required assay in sce
  if (!"normalized.counts" %in% assayNames(sce)) {
    stop("sce must contain a 'normalized.counts' assay.")
  }
  
  
  tryCatch({
    depth.matrix = t(se.var@assays@data$Depth)
    vaf.matrix =   t(se.var@assays@data$VAF)
    variant.ids = rowData(se.var)
    genoqual.matrix = t(se.var@assays@data$Genoqual)
    genotype.matrix = t(se.var@assays@data$Genotype)
    read.counts.df.norm = t(sce@assays@data$normalized.counts)
    cells <- se.var@colData
    if (shiny){
      
      withProgress(message = 'Filter variants', value = 0, {
        incProgress(0, detail = "Depth Filtering...")
        
        # Initial filtering flags
        dp_tf <- depth.matrix < depth.threshold
        incProgress(1/5, detail = "Genotype Quality Filtering...")
        gq_tf <- genoqual.matrix < genotype.quality.threshold
        
        incProgress(1/5, detail = "Variant Allele Frequency and Numerical Genotype Filtering...")
        vaf_ref_tf <- (vaf.matrix > vaf.ref) & (genotype.matrix == 0)
        vaf_hom_tf <- (vaf.matrix < vaf.hom) & (genotype.matrix == 2)
        vaf_het_tf <- (vaf.matrix < vaf.het) & (genotype.matrix == 1)
        
        incProgress(1/5, detail = "Cell and Variant Filtering...")
        # First pass: Determine which entries to keep
        keep <- !(dp_tf | gq_tf | vaf_ref_tf | vaf_hom_tf | vaf_het_tf)
        genotype.matrix[!keep] <- 3
        vaf.matrix[genotype.matrix == 3] <- -1
        
        # Filter based on cell and mutation counts
        num_cells <- nrow(genotype.matrix)
        num_variants <- ncol(genotype.matrix)
        
        cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) {x %in% 0:2})) > num_cells * min.cell / 100
        mut_cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) { x %in% 1:2 })) > num_cells * min.mut.cell / 100
        variant_keep_tf <- cell_num_keep_tf & mut_cell_num_keep_tf  
        v_names <- variant.ids[[1]]
        
        # Second pass filtering
        incProgress(1/5, detail = "Processing Filtered Cells and Variants")
        
        #v_names <- variant.ids
        filtered_variant_names <- v_names[variant_keep_tf]
        cell_variants_keep_tf <- rowSums(genotype.matrix != 3) > num_variants * min.cell / 100
        vaf.matrix.filtered <- vaf.matrix[cell_variants_keep_tf, variant_keep_tf]  
        
        genotype_matrix_filtered <- genotype.matrix[cell_variants_keep_tf, variant_keep_tf]  
        genoqual_matrix_filtered <- genoqual.matrix[cell_variants_keep_tf, variant_keep_tf]
        read.counts.df.norm <- read.counts.df.norm[cell_variants_keep_tf,]
        cells.keep <- cells[cell_variants_keep_tf,]
        
        variant.ids.filtered <- (v_names[variant_keep_tf]) ## Achtung NA
        variant.ids.filtered <- variant.ids.filtered[!is.na(variant.ids.filtered)]
        rownames(vaf.matrix.filtered) <- NULL})
    } else {
      # Initial filtering flags
      dp_tf <- depth.matrix < depth.threshold
      gq_tf <- genoqual.matrix < genotype.quality.threshold
      vaf_ref_tf <- (vaf.matrix > vaf.ref) & (genotype.matrix == 0)
      vaf_hom_tf <- (vaf.matrix < vaf.hom) & (genotype.matrix == 2)
      vaf_het_tf <- (vaf.matrix < vaf.het) & (genotype.matrix == 1)
      
      # First pass: Determine which entries to keep
      keep <- !(dp_tf | gq_tf | vaf_ref_tf | vaf_hom_tf | vaf_het_tf)
      genotype.matrix[!keep] <- 3
      vaf.matrix[genotype.matrix == 3] <- -1
      
      # Filter based on cell and mutation counts
      num_cells <- nrow(genotype.matrix)
      num_variants <- ncol(genotype.matrix)
      cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) {x %in% 0:2})) > num_cells * min.cell / 100
      mut_cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) { x %in% 1:2 })) > num_cells * min.mut.cell / 100
      variant_keep_tf <- cell_num_keep_tf & mut_cell_num_keep_tf  
      v_names <- variant.ids[,1]
      
      # Second pass filtering
      filtered_variant_names <- v_names[variant_keep_tf]
      cell_variants_keep_tf <- rowSums(genotype.matrix != 3) >num_variants  * min.cell / 100
      vaf.matrix.filtered <- vaf.matrix[cell_variants_keep_tf, variant_keep_tf]  
      genotype_matrix_filtered <- genotype.matrix[cell_variants_keep_tf, variant_keep_tf]  
      genoqual_matrix_filtered <- genoqual.matrix[cell_variants_keep_tf, variant_keep_tf]
      #read.counts.df.norm <- read.counts.df.norm[cell_variants_keep_tf,]
      cells.keep <- cells[cell_variants_keep_tf,]
      
      variant.ids.filtered <- (v_names[variant_keep_tf])
      variant.ids.filtered <- variant.ids.filtered[!is.na(variant.ids.filtered)]
      rownames(vaf.matrix.filtered) <- NULL
    }
  }, error = function(e) {
    stop("An error occurred during variant filtering: ", e$message)
  })
  return(list(vaf.matrix.filtered = vaf.matrix.filtered,
              genotype.matrix.filtered = genotype_matrix_filtered,
              #read.counts.df.norm.filtered = read.counts.df.norm,
              variant.ids.filtered = variant.ids.filtered,
              genoqual.matrix.filtered = genoqual_matrix_filtered,
              cells.keep = cells.keep))
}
