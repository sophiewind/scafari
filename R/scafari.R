#' @title scafari: analyzing scDNA-seq data
#'
#' @description scafari is an R Bioconductor package for single-cell DNA-seq 
#' (scDNA-seq) analysis. 
#'
#' @details
#' scafari works on .h5 files, the standard output of the Tapestri pipeline.
#' It offers easy-to-use data quality control as well as explorative variant 
#' analyses and visualization for scDNA-seq data.
#'
#' @section Main Functions:
#' \itemize{
#'   \item{\code{h5ToSce()}: Reads in .h5 files and writes them into a 
#'   SingleCellExperiment object.}
#'   \item{\code{normalizeReadCounts()}: Normalizes the read counts.}
#'   \item{\code{annotateAmplicons()}: Annotates the amplicons  present in the 
#'   .h5 file.}
#'   \item{\code{filterVariants()}: Filters the variants present in the .h5 
#'   file.}
#'   \item{\code{annotateVariants()}: Annotates the filtered variants.}
#'   \item{\code{clusterVariantSelection()}: Clusteres variants of special 
#'   interest.}
#' }
#' 
#' @section Installation:
#' To install this package, use:
#' \code{BiocManager::install("scafari")}
#'
#' @docType package
#' @name scafari
#' @aliases scafari-package
NULL