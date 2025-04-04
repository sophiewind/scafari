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
    
    incProgress(1/5, detail = "Variant Allel Frequency and Numerical Genotype Filtering...")
    vaf_ref_tf <- (vaf.matrix > vaf.ref) & (genotype.matrix == 0)
    vaf_hom_tf <- (vaf.matrix < vaf.hom) & (genotype.matrix == 2)
    vaf_het_tf <- (vaf.matrix < vaf.het) & (genotype.matrix == 1)
    
    incProgress(1/5, detail = "Cell and Variant Filtering...")
    # First pass: Determine which entries to keep
    keep <- !(dp_tf | gq_tf | vaf_ref_tf | vaf_hom_tf | vaf_het_tf)
    #genotype.matrix <- genotype_matrix()
    #vaf.matrix <- vaf_matrix()
    genotype.matrix[!keep] <- 3
    vaf.matrix[genotype.matrix == 3] <- -1
    
    # Filter based on cell and mutation counts
    num_cells <- nrow(genotype.matrix)
    num_variants <- ncol(genotype.matrix)
    
    cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) {x %in% 0:2})) > num_cells * min.cell / 100
    mut_cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) { x %in% 1:2 })) > num_cells * min.mut.cell / 100
    variant_keep_tf <- cell_num_keep_tf & mut_cell_num_keep_tf  
    
    # TODO change
    v_names <- variant.ids$X
    
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
    rownames(vaf.matrix.filtered) <- NULL
    
   
   
    
   
  })
  } else {
    # Initial filtering flags
    dp_tf <- depth.matrix < depth.threshold
    gq_tf <- genoqual.matrix < genotype.quality.threshold
    vaf_ref_tf <- (vaf.matrix > vaf.ref) & (genotype.matrix == 0)
    vaf_hom_tf <- (vaf.matrix < vaf.hom) & (genotype.matrix == 2)
    vaf_het_tf <- (vaf.matrix < vaf.het) & (genotype.matrix == 1)
    
    # First pass: Determine which entries to keep
    keep <- !(dp_tf | gq_tf | vaf_ref_tf | vaf_hom_tf | vaf_het_tf)
    #genotype.matrix <- genotype_matrix()
    #vaf.matrix <- vaf_matrix()
    genotype.matrix[!keep] <- 3
    vaf.matrix[genotype.matrix == 3] <- -1
    
    # Filter based on cell and mutation counts
    num_cells <- nrow(genotype.matrix)
    num_variants <- ncol(genotype.matrix)
    
    cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) {x %in% 0:2})) > num_cells * min.cell / 100
    mut_cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) { x %in% 1:2 })) > num_cells * min.mut.cell / 100
    variant_keep_tf <- cell_num_keep_tf & mut_cell_num_keep_tf  
    
    # TODO change
    v_names <- variant.ids$X
    
    # Second pass filtering
    filtered_variant_names <- v_names[variant_keep_tf]
    cell_variants_keep_tf <- rowSums(genotype.matrix != 3) > num_variants * min.cell / 100
    vaf.matrix.filtered <- vaf.matrix[cell_variants_keep_tf, variant_keep_tf]  
    genotype_matrix_filtered <- genotype.matrix[cell_variants_keep_tf, variant_keep_tf]  
    genoqual_matrix_filtered <- genoqual.matrix[cell_variants_keep_tf, variant_keep_tf]
    read.counts.df.norm <- read.counts.df.norm[cell_variants_keep_tf,]
    cells.keep <- cells[cell_variants_keep_tf,]
    
    variant.ids.filtered <- (v_names[variant_keep_tf])
    variant.ids.filtered <- variant.ids.filtered[!is.na(variant.ids.filtered)]
    rownames(vaf.matrix.filtered) <- NULL
  }
  return(list(vaf.matrix.filtered = vaf.matrix.filtered,
              genotype.matrix.filtered = genotype_matrix_filtered,
              read.counts.df.norm.filtered = read.counts.df.norm,
              variant.ids.filtered = variant.ids.filtered,
              genoqual.matrix.filtered = genoqual_matrix_filtered,
              cells.keep = cells.keep))
}
