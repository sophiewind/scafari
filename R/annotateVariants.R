#' Function: annotateVariants
#' -------------------------------
#' This function takes a SingleCellExperiment object as input and performs variant annotation.
#' 
#' @param sce SingleCellExperiment object containing the single-cell data to be annotated.
#' @param shiny A logical flag indicating whether the function is being run in a Shiny application
#' context. Default is FALSE. 
#' 
#' @return The function returns an annotated SingleCellExperiment object.
#' 
#' @examples
#' 
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp()
#' sce <- annotateVariants(sce, shiny = FALSE)
#'}
#'
#' @export
#' 
#' @references https://missionbio.github.io/mosaic/, https://github.com/rachelgriffard/optima
annotateVariants <- function(sce, shiny = FALSE){
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("The input must be a SingleCellExperiment object.")
  }
  
  # Check if altExp and rowData contain necessary data
  # if (!("X" %in% names(rowData(altExp(sce))))) {
  #   stop("The alternate experiment must have 'X' column in rowData for variant IDs.")
  # }
  # 
  variant.ids.filtered <- rowData(altExp(sce))[[1]]
  metadata = metadata(sce)
  
  # Check that genome version is present in metadata
  if (!"genome_version" %in% names(metadata)) {
    stop("The metadata must contain 'genome_version'.")
  }
  
  if (shiny){
    message('shiny')
    withProgress(message = 'Annotate Variants', value = 0, {
      
      if (check_MBAPI() == 'MissionBio' && metadata[['genome_version']] == 'hg19'){
        # Bring variants to MissionBio APIs variant format
        var.mb <- apply(variant.ids.filtered, 1, function(x){str_match(x, "(chr[0-9XY]+):(\\d+):([ATGC]+)/([ATGC]+)") %>%
            { paste(.[2], .[3], .[4], .[5], sep = "-") }})
        variant.ids.filtered.df.anno <- data.frame(Gene = character(), Protein = character(), `Coding impact` = character(),
                                                   Function = character(), DANN = character(), ClinVar = character() , dbsnp = character())
        
        for (var in seq_along(var.mb)) {
          incProgress(1/length(var.mb), detail = paste0('Annotate variant ', var))
          
          # Reach MissioBio API
          url <- paste0("https://api.missionbio.io/annotations/v1/variants?ids=", var.mb[var])
          
          # Parse URL to get variant annotation
          res = httr::GET(url)
          data = jsonlite::fromJSON(rawToChar(res$content))  # TODO vary
          annot <- data$annotations %>%
            unlist() %>%
            t() %>%
            as.data.frame() 
          
          if (sum (startsWith(colnames(annot), 'function.value')) == 0){
            annot$function.value <- '-'
            
          }
          annot <- annot %>%
            mutate('function.value' = tidyr::unite(
              dplyr::select(., starts_with("function.value")),
              col = 'function.value',
              sep = "-",
              remove = FALSE
            ) %>%
              pull('function.value')) %>% 
            dplyr::select(any_of(c('gene.value', 'protein.value', 'protein_coding_impact.value',
                                   'function.value', 'impact.value', 'clinvar.value', "allele_freq.value", 'dbsnp.value'))) %>%
            dplyr::rename(any_of(c(
              Gene = 'gene.value',
              Protein = 'protein.value',
              `Coding.impact` = 'protein_coding_impact.value',
              Function = 'function.value',
              DANN = 'impact.value',
              ClinVar = 'clinvar.value',
              dbsnp = 'dbsnp.value',
              `Allele Freq (gnomAD)` = "allele_freq.value"))) %>%
            as.data.frame()
          #rownames(annot) <-  var.mb[var]
          variant.ids.filtered.df.anno <- dplyr::bind_rows(variant.ids.filtered.df.anno, annot)  #readRDS('.//input/variant_ids_filtered_df_anno.rds')# 
        }
        
        
        # Original variant ids as rownames
        rownames(variant.ids.filtered.df.anno) <- variant.ids.filtered
        variant.ids.filtered.df.anno$id <- variant.ids.filtered
        
      } else if(check_MBAPI() != 'MissionBio' && metadata[['genome_version']] == 'hg19'){
        stop("MissionBio API is not available.")
      } else if(metadata[['genome_version']] == 'hg38'){
        snpmart <- useEnsembl(biomart = "snp", dataset="hsapiens_snp")
        
        anno <- data.frame('refsnp_source' = c(), 'refsnp_id' = c(), 
                           'chr_name' = c(), 'chrom_start' = c(), 
                           'chrom_end' = c(), 'consequence_type_tv' = c(),
                           'clinical_significance' = c(), 'ensembl_gene_name' = c(), 'id' = c())
        
        for (var in variant.ids.filtered[1:3]){  
          message(paste0('annotate ', var))
          var.tmp <- str_match(var, "([0-9XY]+):(\\d+)") %>% 
            { paste(.[2], .[3], .[3], sep = ":") }
          tryCatch({
            anno.tmp <- getBM(
              attributes = c("refsnp_source", 'refsnp_id', 'chr_name', 'chrom_start', 'chrom_end',
                             "consequence_type_tv", "clinical_significance", 'ensembl_gene_name'),
              filters = 'chromosomal_region',
              values = var.tmp,
              mart = snpmart
            )
            anno.tmp$ID <- var
          }, error = function(e) {
            # Return NA for this SNP if an error occurs
            data.frame(refsnp_source = NA, refsnp_ID = NA, chr_name = NA, chrom_start = NA, 
                       chrom_end = NA, consequence_type_tv = NA, clinical_significance = NA, 
                       ensembl_gene_name = NA)
          })
          anno <- rbind(anno, anno.tmp)
        }
        
        lookup <- bitr(unique(anno$ensembl_gene_name), 'ENSEMBL', 'SYMBOL', org.Hs.eg.db)
        colnames(anno)[colnames(anno) == 'ensembl_gene_name'] <- 'ENSEMBL'
        variant.ids.filtered.df.anno <- merge(anno, lookup, all = T) %>%
          
          dplyr::select(any_of(c('ID', 'SYMBOL', 'Position', 'consequence_type_tv',
                                 'clinical_significance',
                                 'refsnp_id')))%>%
          dplyr::rename(any_of(c(
            Gene = 'SYMBOL',
            `chr_name` = 'Chr',
            `Consequence` = 'consequence_type_tv',
            `Clinical significance` = 'clinical_significance',
            dbsnp = 'refsnp_id'))) %>%
          tidyr::separate(ID, into = c('Chr', 'Start', 'Ref', 'Alt'), sep = ':', remove = F) %>% 
          as.data.frame()
        
      } else {
        message('Issue with MissionBio API. Using Biomart')
      }
      
      
    })
    
    
    # Without shiny logic
  } else {
    message('no shiny')
    
      if (check_MBAPI() == 'MissionBio' && metadata[['genome_version']] == 'hg19'){
        # Bring variants to MissionBio APIs variant format
        var.mb <- apply(variant.ids.filtered, 1, function(x){str_match(x, "(chr[0-9XY]+):(\\d+):([ATGC]+)/([ATGC]+)") %>%
            { paste(.[2], .[3], .[4], .[5], sep = "-") }})
        variant.ids.filtered.df.anno <- data.frame(Gene = character(), Protein = character(), `Coding impact` = character(),
                                                   Function = character(), DANN = character(), ClinVar = character() , dbsnp = character())
        
        for (var in seq_along(var.mb)) {

          # Reach MissioBio API
          url <- paste0("https://api.missionbio.io/annotations/v1/variants?ids=", var.mb[var])
          
          # Parse URL to get variant annotation
          res = httr::GET(url)
          data = jsonlite::fromJSON(rawToChar(res$content))  # TODO vary
          annot <- data$annotations %>%
            unlist() %>%
            t() %>%
            as.data.frame() 
          
          if (sum (startsWith(colnames(annot), 'function.value')) == 0){
            annot$function.value <- '-'
            
          }
          annot <- annot %>%
            mutate('function.value' = tidyr::unite(
              dplyr::select(., starts_with("function.value")),
              col = 'function.value',
              sep = "-",
              remove = FALSE
            ) %>%
              pull('function.value')) %>% 
            dplyr::select(any_of(c('gene.value', 'protein.value', 'protein_coding_impact.value',
                                   'function.value', 'impact.value', 'clinvar.value', "allele_freq.value", 'dbsnp.value'))) %>%
            dplyr::rename(any_of(c(
              Gene = 'gene.value',
              Protein = 'protein.value',
              `Coding.impact` = 'protein_coding_impact.value',
              Function = 'function.value',
              DANN = 'impact.value',
              ClinVar = 'clinvar.value',
              dbsnp = 'dbsnp.value',
              `Allele Freq (gnomAD)` = "allele_freq.value"))) %>%
            as.data.frame()
          #rownames(annot) <-  var.mb[var]
          variant.ids.filtered.df.anno <- dplyr::bind_rows(variant.ids.filtered.df.anno, annot)  #readRDS('.//input/variant_ids_filtered_df_anno.rds')# 
        }
        
        
        # Original variant ids as rownames
        rownames(variant.ids.filtered.df.anno) <- variant.ids.filtered
        variant.ids.filtered.df.anno$id <- variant.ids.filtered
      } else if(check_MBAPI() != 'MissionBio' && metadata[['genome_version']] == 'hg19'){
        # TODO
      } else if(metadata[['genome_version']] == 'hg38'){
        snpmart <- useEnsembl(biomart = "snp", dataset="hsapiens_snp")
        
        anno <- data.frame('refsnp_source' = c(), 'refsnp_id' = c(), 
                           'chr_name' = c(), 'chrom_start' = c(), 
                           'chrom_end' = c(), 'consequence_type_tv' = c(),
                           'clinical_significance' = c(), 'ensembl_gene_name' = c(), 'id' = c())
        
        for (var in variant.ids.filtered[1:3]){  
          message(paste0('annotate ', var))
          var.tmp <- str_match(var, "([0-9XY]+):(\\d+)") %>% 
            { paste(.[2], .[3], .[3], sep = ":") }
          tryCatch({
            anno.tmp <- getBM(
              attributes = c("refsnp_source", 'refsnp_id', 'chr_name', 'chrom_start', 'chrom_end',
                             "consequence_type_tv", "clinical_significance", 'ensembl_gene_name'),
              filters = 'chromosomal_region',
              values = var.tmp,
              mart = snpmart
            )
            anno.tmp$ID <- var
          }, error = function(e) {
            # Return NA for this SNP if an error occurs
            data.frame(refsnp_source = NA, refsnp_ID = NA, chr_name = NA, chrom_start = NA, 
                       chrom_end = NA, consequence_type_tv = NA, clinical_significance = NA, 
                       ensembl_gene_name = NA)
          })
          anno <- rbind(anno, anno.tmp)
        }
        
        lookup <- bitr(unique(anno$ensembl_gene_name), 'ENSEMBL', 'SYMBOL', org.Hs.eg.db)
        colnames(anno)[colnames(anno) == 'ensembl_gene_name'] <- 'ENSEMBL'
        variant.ids.filtered.df.anno <- merge(anno, lookup, all = T) %>%
          
          dplyr::select(any_of(c('ID', 'SYMBOL', 'Position', 'consequence_type_tv',
                                 'clinical_significance',
                                 'refsnp_id')))%>%
          dplyr::rename(any_of(c(
            Gene = 'SYMBOL',
            `chr_name` = 'Chr',
            `Consequence` = 'consequence_type_tv',
            `Clinical significance` = 'clinical_significance',
            dbsnp = 'refsnp_id'))) %>%
          tidyr::separate(ID, into = c('Chr', 'Start', 'Ref', 'Alt'), sep = ':', remove = F) %>% 
          as.data.frame()
        
      } else {
        message('Issue with MissionBio API. Using Biomart')
      }
      
      
  
    
  }
  #rowData(altExp(sce)) <-  variant.ids.filtered.df.anno
  return(variant.ids.filtered.df.anno)
}
