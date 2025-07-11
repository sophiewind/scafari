% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotateAmplicons.R
\name{annotateAmplicons}
\alias{annotateAmplicons}
\title{Function: annotateAmplicons
This function takes a SingleCellExperiment object as input and annotates the
 stored amplicons.}
\usage{
annotateAmplicons(sce, known.canon, shiny = FALSE)
}
\arguments{
\item{sce}{SingleCellExperiment object containing the single-cell data.}

\item{known.canon}{Path to jnown canonicals (see vignette)}

\item{shiny}{If TRUE messages are shown}
}
\value{
A dataframe containing annotated amplicons.
}
\description{
Function: annotateAmplicons
This function takes a SingleCellExperiment object as input and annotates the
 stored amplicons.
}
\examples{
# Assume `sce` is a SingleCellExperiment object with a 'counts' assay
sce_filtered <- readRDS(system.file("extdata", "sce_filtered_demo.rds",
    package = "scafari"
))

annotated <- annotateAmplicons(
    sce_filtered,
    system.file("extdata", "UCSC_hg19_knownCanonical_mock.txt",
        package = "scafari"
    )
)
}
