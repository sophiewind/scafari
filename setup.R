install_if_missing <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    BiocManager::install(package, update = FALSE)
  } else {
    cat(paste0('Package ', package, ' was already installed.\n'))
  }
}

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("BiocInstaller")
BiocInstaller::biocLite('BiocStyle')

# List of packages to install
cran_packages <- c(
  "stringi", "shiny", "shinycssloaders", "DT", "dplyr", "waiter", "ggplot2",
  "tibble", "stringr", "reshape2", "shinyjs", "shinyBS",
  "shinycustomloader", "factoextra", "markdown", "SingleCellExperiment"
)

bioc_packages <- c(
  "GenomicRanges", "rhdf5", "ComplexHeatmap", "karyoploteR",
  "biomaRt", "clusterProfiler", "org.Hs.eg.db"
)

# Install packages
for (package in c(cran_packages, bioc_packages)) {
  tryCatch({
    install_if_missing(package)
  }, error = function(e) {
    cat("Failed to install package:", package, "\n")
    cat("Error message:", conditionMessage(e), "\n")
  })
}

cat("All required packages have been installed.\n")
