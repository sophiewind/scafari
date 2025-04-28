# library(testthat)
# library(scafari)
# 
# # test_that("h5ToSce", {
# #   result <- my_function(known_input)
# #   expect_equal(result, expected_output)
# # })
# 
# test_that("h5ToSce handles errors", {
#   expect_error(h5ToSce('non_existing.h5'), "The file does not exist: non_existing.h5")
# })
# 
# 
# test_that("annotateVariants handles wrong input format", {
#   expect_error(annotateVariants(mtcars), "`sce` must be a SingleCellExperiment object.")
# })
# 
# test_that("clusterVariantSelection handles wrong input format", {
#   expect_error(clusterVariantSelection(mtcars, c()), "`sce` must be a SingleCellExperiment object.")
# })
# 
# # TODO
# test_that("annotateVariants handles wrong input format", {
#   expect_error(filterVariants(mtcars), "`sce` must be a SingleCellExperiment object.")
# })
# 
# test_that("annotateVariants handles wrong input format", {
#   expect_error(logLogPlot(mtcars), "The input must be a SingleCellExperiment object.")
# })
# 
# # NormalizeReadCounts
# test_that("annotateVariants handles wrong input format", {
#   expect_error(normalizeReadCounts(mtcars), "The input must be a SingleCellExperiment object.")
# })
# 
# sce <- readRDS('../inst/extdata/sce_small.rds')
# sce.without.counts <- sce
# sce.without.counts@assays@data$counts <- NULL
# test_that("normalizeCounts handles missing counts", {
#   expect_error(normalizeReadCounts(sce.without.counts), "The SingleCellExperiment object must contain a 'counts' assay.")
# })
# 
# 
# test_that("annotateVariants handles wrong input format", {
#   expect_error(plotAmpliconDistribution(mtcars), "The input must be a SingleCellExperiment object.")
# })

# 
# 
# test_that("plotClusterGenotype handles wrong input format", {
#   expect_error(plotClusterGenotype(mtcars), "The input must be a SingleCellExperiment object.")
# })
# 
# test_that("plotClusterVAFMap handles wrong input format", {
#   expect_error(plotClusterVAFMap(mtcars), "The input must be a SingleCellExperiment object.")
# })
# 
# test_that("plotClusterVAF handles wrong input format", {
#   expect_error(plotClusterVAF(mtcars), "The input must be a SingleCellExperiment object.")
# })
# 
# test_that("plotClusterGenotype handles wrong input format", {
#   expect_error(plotElbow(mtcars, c("FLT3:chr13:28610183:A/G")), "The input must be a SingleCellExperiment object.")
# })

