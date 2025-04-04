readH5 <- function(h5_file){
  # TODO checks
  h5_file <- './input/4-cell-lines-AML-multiomics.dna+protein.h5'
  metadata <- unlist(h5read(h5_file,"assays/dna_read_counts/metadata/"))
  amplicons <- h5read(h5_file, "assays/dna_read_counts/ca/id")
  read.counts <- t(h5read(h5_file, "assays/dna_read_counts/layers/read_counts"))
  
  
  gene.anno.df <- data.frame(
    seqnames = h5read(h5_file, "/assays/dna_read_counts/ca/CHROM"),
    start = h5read(h5_file, "/assays/dna_read_counts/ca/start_pos"),
    end = h5read(h5_file, "/assays/dna_read_counts/ca/end_pos"),
    id = h5read(h5_file, "/assays/dna_read_counts/ca/id")
  )
  if (!startsWith(gene.anno.df$seqnames[1], prefix = 'chr')) gene.anno.df$seqnames <- paste0('chr', gene.anno.df$seqnames)
  
  head(gene.anno.df)
  sce <- SingleCellExperiment(
    assays = list(counts = read.counts),
    colData = gene.anno.df, #amplicons,
    metadata = metadata
  )
}