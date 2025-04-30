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
  
  sample_name <- tryCatch({
    metadata(sce)[['sample_name']]
  }, error = function(e) {
    stop("Failed to extract sample name from metadata: ", e$message)
  })
  
  
  if (genome_version == 'hg19') {
    # Prepare exon database -----------------------------------------------------
    # Read Biomart Exon information and format them
    amps <- as.data.frame(rowData(sce))
    get_exon_data <- function(i) {
      exons <- getBM(
        attributes = c("ensembl_exon_id","ensembl_transcript_id_version", "chromosome_name", "exon_chrom_start", "exon_chrom_end", "rank"),
        filters = c("chromosome_name", "start", "end"),
        values = list(gsub('chr', '', amps$seqnames[i]), amps$start[i], amps$end[i]),
        mart = mart
      )
      exons$region_id <- amps$id[i]  # Include the region ID
      return(exons)
    }
    num_cores <- detectCores() - 1
    
    # Use mclapply for parallel processing
    all_exon_data <- mclapply(1:nrow(amps), get_exon_data, mc.cores = num_cores)
    
    # Combine all exon data into a single data frame
    all_exon_data <- do.call(rbind, all_exon_data)
    head(all_exon_data)
    
    colnames(all_exon_data) <- c('exon_id', 'transcript', 'seqnames', 'start', 'end', 'str', 'rank', 'id')
    all_exon_data$seqnames <- paste0('chr', all_exon_data$seqnames)
    exons.gr <- makeGRangesFromDataFrame(all_exon_data, keep.extra.columns = T)
    
    # Extract canonical transcripts
    canon.path <- system.file("extdata", "UCSC_hg19_knownCanonical_goldenPath.txt", package = "scafari")
    known.canon <- read.delim(canon.path, header = F, col.names = c('seqnames', 'start', 'end', 'x', 'transcript'))
    known.canon$transcript <- gsub('\\..*', '', known.canon$transcript)
    exons.gr$transcript <- gsub('\\..*', '', exons.gr$transcript)
    exons.gr.clean <- exons.gr[exons.gr$transcript %in% known.canon$transcript,]
    gene.anno.gr <- makeGRangesFromDataFrame(amps)
    
    ov <- findOverlaps(gene.anno.gr, exons.gr.clean)
    
    mcols(gene.anno.gr)['transcript'] <- '-'
    gene.anno.gr[queryHits(ov)]$Exon <-     exons.gr.clean[subjectHits(ov)]$rank
    gene.anno.gr[queryHits(ov)]$transcript <-     exons.gr.clean[subjectHits(ov)]$transcript
    # but what if there are multiple????
    # result <- exons.gr.clean %>% 
    #   as.data.frame() %>%
    #   group_by(id) %>%
    #   summarize(
    #     transcripts = paste(unique(transcript), collapse = ", "),
    #     ranks = paste(unique(rank), collapse = ", "),
    #     width = width,
    #     .groups = "drop"  # Ungroups after summarizing
    #   )
    
    # Merge with amps to find unannotated amps
    # amps.anno <- merge(amps, result, by = 'id', all=TRUE) %>% 
    #   replace(is.na(.), '-')
    
    #amps.anno <- amps.anno %>%  dplyr::mutate(Gene = str_split_i(id, '_', 3))
    df <- gene.anno.gr %>%
      as.data.frame()  %>% 
      tibble::rownames_to_column('id') %>% 
      dplyr::mutate(Gene = str_split_i(id, '_', 3)) %>% 
      dplyr::select(seqnames, start, end, width, Gene, Exon, transcript) %>%
      `colnames<-`(c('Chromosome', 'Start', 'End', 'Amplicon length (bp)', 'Gene', 'Exon', 'Canonical Transcript ID')) %>%
      datatable(., rownames = F,  extensions = 'Buttons',
                options = list(pageLength = 10, width = '100%',
                               dom = 'Bfrtip', 
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_panel_", sample_name)),
                                 list(extend = 'excel', filename =  paste0("scafari_panel_",sample_name)),
                                 list(extend = 'pdf', filename =  paste0("scafari_panel_",sample_name)),
                                 list(extend = 'copy', filename =  paste0("scafari_panel_",sample_name)))))
    return(df)
  } else if (genome_version == 'hg38'){
    # MANE annotation
    message('hg38')
    message('MANE annotation is starting. This may take a while.\n')
    mane <- txdbmaker::makeTxDbFromGFF("https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.0/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz")
    mane.df <- exonsBy(mane, by = "tx", use.names = TRUE) %>%  as.data.frame()
    
    message('Processing MANE\n')
    # mane <- mane.gr %>% 
    #   mutate(Exon = str_extract(exon_name, '(?<=exon_number=)\\d+'), 
    #         # Gene = str_extract(V9, '(?<=gene_name=)[^;]+'),  
    #         # `Transcript ID` = str_extract(V9, '(?<=transcript_id=)[^;]+')
    #         )  %>% 
    #   filter(V3 == 'exon')
    # colnames(mane) <- c('seqnames','source', 'feature','start', 'end','score','strand', 'frame','Atrribute', 'Exon', 'Gene', 'Transcript ID')
    
    message('Annotating\n')
    mane.gr <- makeGRangesFromDataFrame(mane.df, keep.extra.columns = T, na.rm = T)
    
    # Define transcript columns
    gene.anno.gr <- makeGRangesFromDataFrame(rowData(sce))
    mcols(gene.anno.gr)[["Exon"]] <- '-'
    mcols(gene.anno.gr)[["Transcript ID"]] <- '-'
    mcols(gene.anno.gr)[["Gene"]] <- '-'
    ov <- findOverlaps(gene.anno.gr, mane.gr)
    
    # Until here done
    mcols(gene.anno.gr)[queryHits(ov),][["Exon"]] <- mcols(mane.gr)[subjectHits(ov),][["exon_rank"]]
    mcols(gene.anno.gr)[queryHits(ov),]["Transcript ID"] <- mcols(mane.gr)[subjectHits(ov),]["group_name"]
    #mcols(gene.anno.gr)[queryHits(ov),]["Gene"] <- mcols(mane.gr)[subjectHits(ov),]["Gene"]
    
    df <- gene.anno.gr %>% as.data.frame() %>% 
      dplyr::select(id, seqnames, start, end, width, Gene, Exon, `Transcript.ID`) %>% 
      `colnames<-`(c('Amplicon ID', 'Chromosome', 'Start', 'End', 'Amplicon length (bp)', 'Gene', 'Exon', 'Canonical Transcript ID')) %>% 
      datatable(.,  extensions = 'Buttons',
                options = list(pageLength = 10, width = '100%',
                               dom = 'Bfrtip', 
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_panel_",sample_name)),
                                 list(extend = 'excel', filename =  paste0("scafari_panel_",sample_name)),
                                 list(extend = 'pdf', filename =  paste0("scafari_panel_",sample_name)),
                                 list(extend = 'copy', filename =  paste0("scafari_panel_",sample_name)))), rownames = F)
    return(df)
  } else {
    message('No proper genome version')
  }
}
