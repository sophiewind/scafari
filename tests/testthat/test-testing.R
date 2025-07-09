# Setup data -------------------------------------------------------------------
# Create test dataset to manipulate
library(SingleCellExperiment)
h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
h5 <- h5ToSce(h5_file_path)
sce <- h5$sce_amp
se.var <- h5$se_var

sce_filtered <- readRDS(system.file("extdata", "sce_filtered_demo.rds",
    package = "scafari"
))
sce_filtered.wrong_genome <- sce_filtered
metadata(sce_filtered.wrong_genome)[["genome_version"]] <- "test"

sce_filtered.without_vaf <- sce_filtered
SingleCellExperiment::altExp(sce_filtered.without_vaf)[["VAF"]] <- NULL

clusterplot <- readRDS(system.file("extdata", "clusterplot.rds",
    package = "scafari"
))
# Remove normalized counts
sce.without.counts <- sce
sce.without.counts@assays@data$counts <- NULL

# Remove normalized counts
sce.without.meta <- sce
sce.without.meta@metadata <- list(c())

# remove VAF assay
se.var.missing <- se.var
assays(se.var.missing)[["VAF"]] <- NULL

# Error testing ----------------------------------------------------------------
## IO --------------------------------------------------------------------------
test_that("h5ToSce handles errors", {
    expect_error(
        h5ToSce("non_existing.h5"),
        "The file does not exist: non_existing.h5"
    )
})

## Sequencing analysis ---------------------------------------------------------
test_that("logLogPlot handles wrong input format", {
    expect_error(
        logLogPlot(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("logLogPlot handles wrong input format", {
    expect_error(
        logLogPlot(sce.without.counts),
        "The SingleCellExperiment object must contain a 'counts' assay."
    )
})

## Panel analysis --------------------------------------------------------------
test_that("normalizeReadCounts handles wrong input format", {
    expect_error(
        normalizeReadCounts(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("normalizeCounts handles missing counts", {
    expect_error(
        normalizeReadCounts(sce.without.counts),
        "The SingleCellExperiment object must contain a 'counts' assay."
    )
})

test_that("normalizeCounts handles missing metadata", {
    expect_error(
        normalizeReadCounts(sce.without.meta),
        "The metadata must contain 'n_amplicons'."
    )
})

test_that("plotNormalizedReadCounts handles wrong input format", {
    expect_error(
        plotNormalizedReadCounts(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("plotNormalizedReadCounts handles missing norm.counts", {
    expect_error(
        plotNormalizedReadCounts(sce.without.counts),
        "The SingleCellExperiment object must contain a 'normalized.counts' assay."
    )
})

test_that("plotAmpliconDistribution handles wrong input format", {
    expect_error(
        plotAmpliconDistribution(sce = mtcars),
        "The input must be a SingleCellExperiment object."
    )
})


test_that("plotPanelUniformity handles wrong input format", {
    expect_error(
        plotPanelUniformity(mtcars),
        "`sce` must be a SingleCellExperiment object."
    )
})

test_that("plotPanelUniformity handles wrong input format", {
    expect_error(
        plotPanelUniformity(sce, interactive = "test"),
        "`interactive` must be logical."
    )
})

test_that("plotPanelUniformity handles wrong input format", {
    expect_error(
        plotPanelUniformity(sce.without.counts),
        "The SingleCellExperiment object must contain a 'normalized.counts' assay."
    )
})


test_that("annotateAmplicons handles wrong input format", {
    expect_error(
        annotateAmplicons(sce = mtcars),
        "`sce` must be a SingleCellExperiment object."
    )
})

test_that("annotateAmplicons handles missing genome_version", {
    expect_error(
        annotateAmplicons(sce.without.meta),
        "The metadata must contain 'genome_version'."
    )
})

test_that("annotateAmplicons handles missing genome_version", {
    expect_error(
        annotateAmplicons(sce.without.meta),
        "The metadata must contain 'genome_version'."
    )
})

test_that("plotGenotypequalityPerGenotype handles invalid inputs properly", {
    expect_error(
        plotGenotypequalityPerGenotype(sce = mtcars),
        "Input 'sce' is invalid. It must be a SingleCellExperiment object."
    )
})

test_that("annotateAmplicons handles missing metadata", {
    expect_error(
        annotateAmplicons(sce.without.meta),
        "The metadata must contain 'genome_version'."
    )
})


## Variant analysis ------------------------------------------------------------
test_that("filterVariants handles wrong input format", {
    expect_error(
        filterVariants(
            depth.threshold = 10,
            genotype.quality.threshold = 30,
            vaf.ref = 5,
            vaf.het = 35,
            vaf.hom = 95,
            min.cell = 50,
            min.mut.cell = 1,
            se.var = mtcars,
            sce = sce,
            shiny = FALSE
        ),
        "se.var must be a SummarizedExperiment object."
    )
})


test_that("filterVariants handles wrong input format", {
    expect_error(
        filterVariants(
            depth.threshold = 10,
            genotype.quality.threshold = 30,
            vaf.ref = 5,
            vaf.het = 35,
            vaf.hom = 95,
            min.cell = 50,
            min.mut.cell = 1,
            se.var = se.var,
            sce = se.var,
            shiny = FALSE
        ),
        "sce must be a SingleCellExperiment object."
    )
})

test_that("filterVariants handles wrong input format", {
    expect_error(
        filterVariants(
            depth.threshold = 10,
            genotype.quality.threshold = 30,
            vaf.ref = 5,
            vaf.het = 35,
            vaf.hom = 95,
            min.cell = 50,
            min.mut.cell = 1,
            se.var = se.var.missing,
            sce = sce,
            shiny = FALSE
        ),
        "Missing required assays in se.var: VAF"
    )
})



test_that("filterVariants handles wrong input format", {
    expect_error(
        filterVariants(
            depth.threshold = 10,
            genotype.quality.threshold = 30,
            vaf.ref = 5,
            vaf.het = 35,
            vaf.hom = 95,
            min.cell = 50,
            min.mut.cell = 1,
            se.var = se.var,
            sce = sce.without.counts,
            shiny = FALSE
        ),
        "sce must contain a 'normalized.counts' assay."
    )
})


test_that("annotateVariants handles wrong input format", {
    expect_error(
        annotateVariants(mtcars),
        "`sce` must be a SingleCellExperiment object."
    )
})


test_that("annotateVariants handles wrong input format", {
    expect_error(
        annotateVariants(sce, max.var = "test"),
        "`max.var` must be a numeric."
    )
})

test_that("annotateVariants handles wrong input format", {
    expect_error(
        annotateVariants(sce, shiny = "test"),
        "`shiny` must be logical."
    )
})

test_that("annotateVariants handles wrong max.var format", {
    expect_error(
        annotateVariants(sce, shiny = FALSE, max.var = "test"),
        "`max.var` must be a numeric."
    )
})


test_that("plotVariantHeatmap handles wrong input format", {
    expect_error(
        plotVariantHeatmap(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("plotVariantHeatmap handles missing altExp", {
    expect_error(
        plotVariantHeatmap(sce),
        "The SingleCellExperiment object must contain 'variants' as an alternate experiment."
    )
})


## Explore variants analysis ---------------------------------------------------
test_that("plotElbow handles wrong input format", {
    expect_error(
        plotElbow(sce = mtcars, variants.of.interest = c()),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("plotElbow handles wrong input format", {
    expect_error(
        plotElbow(sce = sce, variants.of.interest = c()),
        "The 'variants.of.interest' parameter must be a vector."
    )
})

test_that("clusterVariantSelection handles wrong input format", {
    expect_error(
        clusterVariantSelection(mtcars, c()),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("clusterVariantSelection handles wrong input format", {
    expect_error(
        clusterVariantSelection(sce, c()),
        "The SingleCellExperiment object must contain 'variants' as an alternate experiment."
    )
})

test_that("clusterVariantSelection handles wrong input variarbles", {
    expect_error(
        clusterVariantSelection(sce_filtered, c("test")),
        "All variants.of.interest must exist in the VAF matrix columns."
    )
})

test_that("plotClusterGenotype handles wrong input format", {
    expect_error(
        plotClusterGenotype(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("plotClusterGenotype handles empty variants vector.", {
    expect_error(
        plotClusterGenotype(sce_filtered, c()),
        "variants.of.interest must be a non-empty vector."
    )
})

test_that("plotClusterGenotype handles wrong variables", {
    expect_error(
        plotClusterGenotype(sce_filtered, c("test_variant")),
        "All variants.of.interest must exist in the VAF matrix columns."
    )
})

test_that("plotClusterVAFMap handles wrong input format", {
    expect_error(
        plotClusterVAFMap(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("plotClusterVAFMap handles empty variants vector.", {
    expect_error(
        plotClusterVAFMap(sce_filtered, c()),
        "variants.of.interest must be a non-empty vector."
    )
})

test_that("plotClusterGenotype handles wrong variables", {
    expect_error(
        plotClusterVAFMap(
            sce_filtered, c("test_variant"),
            clusterplot$clusterplot
        ),
        "All variants.of.interest must exist in the VAF matrix columns."
    )
})


test_that("plotClusterVAF handles wrong input format", {
    expect_error(
        plotClusterVAF(mtcars),
        "The input must be a SingleCellExperiment object."
    )
})

test_that("plotClusterVAF handles empty variants vector.", {
    expect_error(
        plotClusterVAF(sce_filtered, c()),
        "variants.of.interest must be a non-empty vector."
    )
})

test_that("plotClusterVAF handles wrong variables", {
    expect_error(
        plotClusterVAF(sce_filtered, c("test_variant")),
        "All variants.of.interest must exist in the VAF matrix columns."
    )
})


test_that("plotClusterGenotype handles wrong input format", {
    expect_error(
        plotElbow(mtcars, c("FLT3:chr13:28610183:A/G")),
        "The input must be a SingleCellExperiment object."
    )
})

# Function tests ---------------------------------------------------------------
test_that("h5ToSce creates correct output.", {
    expect_equal(h5ToSce(h5_file_path), h5)
})

# sce
test_that("normalizeReadCounts Normalization of counts is working.", {
    norm <- normalizeReadCounts(sce)
    sce_test <- readRDS(system.file("extdata", "sce_norm_demo.rds",
        package = "scafari"
    ))
    expect_equal(norm@assays@data$normalized.counts, sce_test@assays@data$normalized.counts)
})

test_that("annotateAmplicons", {
    # Since its only possible to annotate by mock in test not all columns are
    # considered
    annotated <- annotateAmplicons(
        sce,
        system.file("extdata",
            "UCSC_hg19_knownCanonical_mock.txt",
            package = "scafari"
        )
    )
    annotated.true <- readRDS(system.file("extdata",
        "annotated.rds",
        package = "scafari"
    ))
    expect_equal(annotated$x$data[3:5, ], annotated.true$x$data[3:5, ])
})


test_that("filterVariants creates correct output.", {
    h5 <- h5ToSce(h5_file_path)
    sce <- h5$sce_amp
    se.var <- h5$se_var
    sce <- normalizeReadCounts(sce)
    filteres <- filterVariants(
        depth.threshold = 10,
        genotype.quality.threshold = 30,
        vaf.ref = 5,
        vaf.het = 35,
        vaf.hom = 95,
        min.cell = 50,
        min.mut.cell = 1,
        se.var = se.var,
        sce = sce,
        shiny = FALSE
    )
    se.f <-
        SummarizedExperiment(
            assays =
                list(
                    VAF = t(filteres$vaf.matrix.filtered),
                    Genotype = t(filteres$genotype.matrix.filtered),
                    Genoqual =
                        t(filteres$genoqual.matrix.filtered)
                ),
            rowData = filteres$variant.ids.filtered,
            colData = filteres$cells.keep
        )

    # Filter out cells in sce object
    # Find the indices of the columns to keep
    indices_to_keep <- match(filteres$cells.keep,
        SummarizedExperiment::colData(sce)[[1]],
        nomatch = 0
    )

    # Subset the SCE using these indices
    sce_filtered_test <- sce[, indices_to_keep]
    SingleCellExperiment::altExp(sce_filtered_test, "variants") <- se.f
    # sce_filtered_test <- annotateVariants(sce = sce_filtered)

    # Because APIs can differ its only possible to check parts of the output
    expect_equal(
        assays(SingleCellExperiment::altExp(sce_filtered_test, "variants"))[["VAF"]],
        assays(SingleCellExperiment::altExp(sce_filtered, "variants"))[["VAF"]]
    )
})


# sce
test_that("Normalization of counts is working.", {
    norm <- normalizeReadCounts(sce)
    sce_test <- readRDS(system.file("extdata", "sce_norm_demo.rds",
        package = "scafari"
    ))
    expect_equal(norm, sce_test)
})

test_that("filterVariants creates correct output.", {
    h5 <- h5ToSce(h5_file_path)
    sce <- h5$sce_amp
    se.var <- h5$se_var
    sce <- normalizeReadCounts(sce)
    filteres <- filterVariants(
        depth.threshold = 10,
        genotype.quality.threshold = 30,
        vaf.ref = 5,
        vaf.het = 35,
        vaf.hom = 95,
        min.cell = 50,
        min.mut.cell = 1,
        se.var = se.var,
        sce = sce,
        shiny = FALSE
    )
    se.f <-
        SummarizedExperiment(
            assays =
                list(
                    VAF = t(filteres$vaf.matrix.filtered),
                    Genotype = t(filteres$genotype.matrix.filtered),
                    Genoqual =
                        t(filteres$genoqual.matrix.filtered)
                ),
            rowData = filteres$variant.ids.filtered,
            colData = filteres$cells.keep
        )

    # Filter out cells in sce object
    # Find the indices of the columns to keep
    indices_to_keep <- match(filteres$cells.keep,
        SummarizedExperiment::colData(sce)[[1]],
        nomatch = 0
    )

    # Subset the SCE using these indices
    sce_filtered_test <- sce[, indices_to_keep]
    SingleCellExperiment::altExp(sce_filtered_test, "variants") <- se.f
    # sce_filtered_test <- annotateVariants(sce = sce_filtered)

    # Because APIs can differ its only possible to check parts of the output
    expect_equal(
        assays(SingleCellExperiment::altExp(sce_filtered_test, "variants"))[["VAF"]],
        assays(SingleCellExperiment::altExp(sce_filtered, "variants"))[["VAF"]]
    )
})


test_that("Variant annotation works", {
    # Because different APIs can lead to different columns, just the presence of
    # the id column can be tested
    annotated.var <- annotateVariants(sce_filtered)
    expect_true("id" %in% colnames(rowData(
        SingleCellExperiment::altExp(sce_filtered)
    )))
})

