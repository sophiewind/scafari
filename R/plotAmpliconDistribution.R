plotAmpliconDistribution <- function(sce) {
  gene_anno_df <- as.data.frame(rowData(sce))
  # Ensure the input data is a data frame
  if (!is.data.frame(gene_anno_df)) {
    stop("gene_anno_df must be a data frame.")
  }
  
  # Prepare GRanges from gene annotation data frame
  gene_anno_gr <- gene_anno_df %>%
    mutate(Gene = sapply(strsplit(id, '_'), function(x) x[3])) %>% # Assuming str_split_i() splits by '_'
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # Create a karyotype plot
  kp <- plotKaryotype(plot.type = 1, cex = 1.5)
  
  # Plot genomic regions
  kpPlotRegions(kp, data = gene_anno_gr, col = "#F66D7A", data.panel = 1)
  
  # Prepare GRanges object for gene labels (avoid overplotting)
  amp_genes <- gene_anno_df %>%
    mutate(Gene = sapply(strsplit(id, '_'), function(x) x[3])) %>%
    group_by(Gene) %>%
    slice(1) %>%
    ungroup() %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)
  
  # Add gene labels to the plot
  kpText(kp, chr = seqnames(amp_genes), x = start(amp_genes), y = 0.75, labels = amp_genes$Gene, cex = 0.95)
}
