# Setup data -------------------------------------------------------------------
# Create test dataset to manipulate
h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
h5 <- h5ToSce(h5_file_path)
sce <- h5$sce_amp
se.var <- h5$se_var

# Remove normalized counts
sce.without.counts <- sce
sce.without.counts@assays@data$counts <- NULL

# remove VAF assay
se.var.missing <- se.var
assays(se.var.missing)[['VAF']] <- NULL

# IO ---------------------------------------------------------------------------
test_that("h5ToSce handles errors", {
  expect_error(h5ToSce('non_existing.h5'), 
               "The file does not exist: non_existing.h5")
})

# Panel analysis ---------------------------------------------------------------
test_that("annotateVariants handles wrong input format", {
  expect_error(logLogPlot(mtcars), 
               "The input must be a SingleCellExperiment object.")
})

test_that("annotateVariants handles wrong input format", {
  expect_error(normalizeReadCounts(mtcars), 
               "The input must be a SingleCellExperiment object.")
})


test_that("normalizeCounts handles missing counts", {
  expect_error(normalizeReadCounts(sce.without.counts), 
               "The SingleCellExperiment object must contain a 'counts' assay.")
})

test_that("annotateVariants handles wrong input format", {
  expect_error(plotAmpliconDistribution(mtcars), 
               "The input must be a SingleCellExperiment object.")
})


# Variant analysis -------------------------------------------------------------
test_that("annotateVariants handles wrong input format", {
  expect_error(filterVariants(depth.threshold = 10,
                              genotype.quality.threshold = 30,
                              vaf.ref = 5, 
                              vaf.het =  35, 
                              vaf.hom = 95, 
                              min.cell = 50,
                              min.mut.cell = 1,
                              se.var = mtcars,
                              sce = sce,
                              shiny = FALSE), 
               "se.var must be a SummarizedExperiment object.")
})


test_that("annotateVariants handles wrong input format", {
  expect_error(filterVariants(depth.threshold = 10,
                              genotype.quality.threshold = 30,
                              vaf.ref = 5, 
                              vaf.het =  35, 
                              vaf.hom = 95, 
                              min.cell = 50,
                              min.mut.cell = 1,
                              se.var = se.var,
                              sce = se.var,
                              shiny = FALSE), 
               "sce must be a SingleCellExperiment object.")
})

test_that("annotateVariants handles wrong input format", {
  expect_error(filterVariants(depth.threshold = 10,
                              genotype.quality.threshold = 30,
                              vaf.ref = 5, 
                              vaf.het =  35, 
                              vaf.hom = 95, 
                              min.cell = 50,
                              min.mut.cell = 1,
                              se.var = se.var.missing,
                              sce = sce,
                              shiny = FALSE), 
               "Missing required assays in se.var: VAF")
})



test_that("annotateVariants handles wrong input format", {
  expect_error(filterVariants(depth.threshold = 10,
                              genotype.quality.threshold = 30,
                              vaf.ref = 5, 
                              vaf.het =  35, 
                              vaf.hom = 95, 
                              min.cell = 50,
                              min.mut.cell = 1,
                              se.var = se.var,
                              sce = sce.without.counts,
                              shiny = FALSE), 
               "sce must contain a 'normalized.counts' assay.")
})


test_that("annotateVariants handles wrong input format", {
  expect_error(annotateVariants(mtcars), 
               "`sce` must be a SingleCellExperiment object.")
})

# Explore variants analysis ----------------------------------------------------
test_that("clusterVariantSelection handles wrong input format", {
  expect_error(clusterVariantSelection(mtcars, c()), 
               "The input must be a SingleCellExperiment object.")
})


test_that("plotClusterGenotype handles wrong input format", {
  expect_error(plotClusterGenotype(mtcars), 
               "The input must be a SingleCellExperiment object.")
})

test_that("plotClusterVAFMap handles wrong input format", {
  expect_error(plotClusterVAFMap(mtcars), 
               "The input must be a SingleCellExperiment object.")
})

test_that("plotClusterVAF handles wrong input format", {
  expect_error(plotClusterVAF(mtcars), 
               "The input must be a SingleCellExperiment object.")
})

test_that("plotClusterGenotype handles wrong input format", {
  expect_error(plotElbow(mtcars, c("FLT3:chr13:28610183:A/G")), 
               "The input must be a SingleCellExperiment object.")
})

