#' Function: annotateAmplicons
#' -------------------------------
#' This function takes a SingleCellExperiment object as input and annotates the stored amplicons.
#' 
#' @param sce SingleCellExperiment object containing the single-cell data.
#' 
#' @return A dataframe containing annotated amplicons.
#' 
#' \dontrun{
#' # Assume `sce` is a SingleCellExperiment object with a 'counts' assay
#' annotated <- annotateAmplicons(sce)
#'}
#'
#' @export
annotateAmplicons <- function(sce){
  # Check if the SCE object has metadata
  if (is.null(metadata(sce))) {
    stop("The SingleCellExperiment object does not contain any metadata.")
  }
  
  # Attempt to extract genome version
  genome_version <- tryCatch({
    metadata(sce)[['genome_version']]
  }, error = function(e) {
    stop("Failed to extract genome version from metadata: ", e$message)
  })
  
  if (genome_version == 'hg19') {
    # Prepare exon database -----------------------------------------------------
    # Read Biomart Exon information and format them
    exons_path <- './input/mart_export_grch37_p13.txt'
    
    # Check if the file exists
    if (!file.exists(exons_path)) {
      stop("The exon file does not exist at path: ", exons_path)
    }
    
    exons <- tryCatch({
      read.delim(exons_path, header = TRUE, sep = '\t')
    }, error = function(e) {
      stop("Failed to read the exon information file: ", e$message)
    })
    
    colnames(exons)[11] <- 'seqnames'
    colnames(exons)[5] <- 'start'
    colnames(exons)[6] <- 'end'
    colnames(exons)[7] <- 'Exon'
    exons$seqnames <- paste0('chr', exons$seqnames)
    exons.gr <- makeGRangesFromDataFrame(exons, keep.extra.columns = T)
    
    # Extract canonical transcripts
    known.canon <- read.delim('./input/UCSC_hg19_knownCanonical_chrom.bed', header = F, col.names = c('seqnames', 'start', 'end', 'Transcript.stable.ID.version', 'gene'))
    known.canon$Transcript.stable.ID.version <- gsub('\\..*', '', known.canon$Transcript.stable.ID.version)
    exons.gr$Transcript.stable.ID.version <- gsub('\\..*', '', exons.gr$Transcript.stable.ID.version)
    exons.gr.clean <- exons.gr[exons.gr$Transcript.stable.ID.version %in% known.canon$Transcript.stable.ID.version,]
    
    # Change chr info to merge them
    gene.anno.df.tmp <- gene.anno.df
    #gene.anno.df.tmp$seqnames <- gsub('chr', '', gene.anno.df()$seqnames)
    gene.anno.df.tmp <- gene.anno.df.tmp %>% dplyr::mutate(Gene = str_split_i(id, '_', 3))
    gene.anno.gr <- makeGRangesFromDataFrame(gene.anno.df.tmp, keep.extra.columns = T)
    rm(gene.anno.df.tmp)
    exons.gr.clean <- exons.gr[exons.gr$Transcript.stable.ID.version %in%
                                 known.canon$Transcript.stable.ID.version,]
    ov <- findOverlaps(gene.anno.gr, exons.gr.clean, )
    
    # Define transcript columns
    mcols(gene.anno.gr)$`Exon` <- '-'
    mcols(gene.anno.gr)$`Canonical Transcript ID` <- '-'
    
    gene.anno.gr[queryHits(ov)]$`Exon` <- exons.gr.clean[subjectHits(ov)]$Exon  
    gene.anno.gr[queryHits(ov)]$`Canonical Transcript ID` <- exons.gr.clean[subjectHits(ov)]$Transcript.stable.ID.version  # i think its wrong
    
    # Format for datatable
    df <- gene.anno.gr %>% as.data.frame() %>% 
      dplyr::select(seqnames, start, end, width, Gene, Exon, Canonical.Transcript.ID) %>%
      `colnames<-`(c('Chromosome', 'Start', 'End', 'Amplicon length (bp)', 'Gene', 'Exon', 'Canonical Transcript ID')) %>%
      datatable(., rownames = F,  extensions = 'Buttons',
                options = list(pageLength = 10, width = '100%',
                               dom = 'Bfrtip', 
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_panel_", metadata[['sample_name']])),
                                 list(extend = 'excel', filename =  paste0("scafari_panel_", metadata[['sample_name']])),
                                 list(extend = 'pdf', filename =  paste0("scafari_panel_", metadata[['sample_name']])),
                                 list(extend = 'copy', filename =  paste0("scafari_panel_", metadata[['sample_name']])))))
    return(df)
  } else if (metadata[['genome_version']] == 'hg38'){
    # MANE annotation
    message('hg38')
    message('MANE annotation is starting. This may take a while.\n')
    mane.raw <- read.delim('./input/MANE.GRCh38.v1.3.ensembl_genomic.gff', skip = 2, header = F, sep = '\t')
    
    message('Processing MANE\n')
    mane <- mane.raw %>% 
      mutate(Exon = str_extract(V9, '(?<=exon_number=)\\d+'), 
             Gene = str_extract(V9, '(?<=gene_name=)[^;]+'),  
             `Transcript ID` = str_extract(V9, '(?<=transcript_id=)[^;]+'))  %>% 
      filter(V3 == 'exon')
    colnames(mane) <- c('seqnames','source', 'feature','start', 'end','score','strand', 'frame','Atrribute', 'Exon', 'Gene', 'Transcript ID')
    
    message('Annotating\n')
    mane.gr <<- makeGRangesFromDataFrame(mane, keep.extra.columns = T, na.rm = T)
    
    # Define transcript columns
    gene.anno.gr <- gene.anno.gr
    mcols(gene.anno.gr)[["Exon"]] <- '-'
    mcols(gene.anno.gr)[["Transcript ID"]] <- '-'
    mcols(gene.anno.gr)[["Gene"]] <- '-'
    ov <- findOverlaps(gene.anno.gr, mane.gr)
    
    mcols(gene.anno.gr)[queryHits(ov),]["Exon"] <- mcols(mane.gr)[subjectHits(ov),]["Exon"]
    mcols(gene.anno.gr)[queryHits(ov),]["Transcript ID"] <- mcols(mane.gr)[subjectHits(ov),]["Transcript ID"]
    mcols(gene.anno.gr)[queryHits(ov),]["Gene"] <- mcols(mane.gr)[subjectHits(ov),]["Gene"]
    
    df <- gene.anno.gr %>% as.data.frame() %>% 
      dplyr::select(id, seqnames, start, end, width, Gene, Exon, `Transcript.ID`) %>% 
      `colnames<-`(c('Amplicon ID', 'Chromosome', 'Start', 'End', 'Amplicon length (bp)', 'Gene', 'Exon', 'Canonical Transcript ID')) %>% 
      datatable(.,  extensions = 'Buttons',
                options = list(pageLength = 10, width = '100%',
                               dom = 'Bfrtip', 
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_panel_", metadata[['sample_name']])),
                                 list(extend = 'excel', filename =  paste0("scafari_panel_", metadata[['sample_name']])),
                                 list(extend = 'pdf', filename =  paste0("scafari_panel_", metadata[['sample_name']])),
                                 list(extend = 'copy', filename =  paste0("scafari_panel_", metadata[['sample_name']])))), rownames = F)
    return(df)
  } else {
    message('No proper genome version')
  }
}
