#' Function: h5ToSce
#' 
#' This function takes the path of an h5 file and reads this into an object of the SingleCellExperiment class.
#' 
#' @param h5_file  path of an h5 file
#' 
#' @return A list containing the following elements:
#' \describe{
#'   \item{sce_amp}{SingleCellExperiment class object containing read count information.}
#'   \item{se_var}{SummarizedExperiment class object containing variant information..}
#' } 
#' 
#' 
#' @examples
#' h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
#'
#' # Read the h5ToSce using readH5File
#' result <- h5ToSce(h5_file_path)
#'
#' # Display the result
#' print(result)
#' 
#' @return A list with the SingleCellExperiment and the SummarizedExperiment 
#' object representing the amplicon and the varaint analysis.
#' \describe{
#'   \item{sce_amp}{SingleCellExperiment with amplicon experiment.}
#'   \item{se_var}{SummarizedExperiment with variants eexperiment}
#' }
#' 
#' @import SummarizedExperiment
#' 
#' @export
h5ToSce <- function(h5_file) {
  # Process amplicon data ------------------------------------------------------
  # Check if the file exists
  if (!file.exists(h5_file)) {
    stop("The file does not exist: ", h5_file)
  }
  
  # Try reading the metadata, stop if it fails
  tryCatch({
    metadata <- unlist(h5read(h5_file, "assays/dna_variants/metadata/"))
  }, error = function(e) {
    stop("Failed to read metadata: ", e$message)
  })
  
  # Try reading variant IDs, stop if it fails
  tryCatch({
    variant.ids <- h5read(h5_file, 'assays/dna_variants/ca/id')
  }, error = function(e) {
    stop("Failed to read variant IDs: ", e$message)
  })
  

  # Try reading cell barcodes, stop if it fails
  tryCatch({
    cells.rc <- h5read(h5_file, "assays/dna_read_counts/ra/barcode")
    cells.var <- h5read(h5_file, "assays/dna_variants/ra/barcode")
  }, error = function(e) {
    stop("Failed to read cell barcode information: ", e$message)
  })
  
  # Try reading depth matrix and other matrices, stop if it fails
  tryCatch({
    depth.matrix <-h5read(h5_file, "assays/dna_variants/layers/DP")
    genoqual.matrix <- h5read(h5_file, "assays/dna_variants/layers/GQ")
    genotype.matrix <- h5read(h5_file, "assays/dna_variants/layers/NGT")
    vaf.matrix <- h5read(h5_file, "assays/dna_variants/layers/AF")
    amplicons <- h5read(h5_file, 'assays/dna_read_counts/ca/id')
  }, error = function(e) {
    stop("Failed to read variant or associated matrices: ", e$message)
  })
  
  # Try reading read counts, stop if it fails
  tryCatch({
    read.counts.df <- as.data.frame(t(h5read(h5_file, "assays/dna_read_counts/layers/read_counts")))
    colnames(read.counts.df) <- amplicons
  }, error = function(e) {
    stop("Failed to read read counts: ", e$message)
  })
  
  # Try reading gene annotation, stop if it fails
  tryCatch({
    gene.anno.df <- data.frame(
      seqnames = paste0('chr', h5read(h5_file, "/assays/dna_read_counts/ca/CHROM")),
      start = h5read(h5_file, "/assays/dna_read_counts/ca/start_pos"),
      end = h5read(h5_file, "/assays/dna_read_counts/ca/end_pos"),
      id = h5read(h5_file, "/assays/dna_read_counts/ca/id")
    )
  }, error = function(e) {
    stop("Failed to read gene annotations: ", e$message)
  })
  
  # Create the SingleCellExperiment object when data is ready
  tryCatch({
    sce <- SingleCellExperiment(
      assays = list(counts = t(read.counts.df)),
      rowData = gene.anno.df,
      metadata = metadata,
      colData = cells.rc
    )
  }, error = function(e) {
    stop("Failed to create SingleCellExperiment object: ", e$message)
  })
  
  # Process variant data -------------------------------------------------------
  se <- SummarizedExperiment(assays = list(VAF = vaf.matrix, 
                                     Genotype = genotype.matrix,
                                     Genoqual = genoqual.matrix,
                                     Depth = depth.matrix
                                     ),
                       rowData  = DataFrame(variant.ids),
                       colData = cells.var)
  return(list(sce_amp = sce, se_var = se))
}
