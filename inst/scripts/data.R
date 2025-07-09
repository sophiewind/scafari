# Documentation how inst/extdata was generated

# annotated.rds ----------------------------------------------------------------
# h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
# h5 <- h5ToSce(h5_file_path)
# sce <- h5$sce_amp
# se.var <- h5$se_var
# sce <- normalizeReadCounts(sce = sce)
# annotated <- try(annotateAmplicons(sce = sce, known.canon = known_canon_path))
# saveRDS(annotated, './inst/extdata/annotated.rds')

# demo.h5 ----------------------------------------------------------------------
# .h5 was downloaded here: 
# https://portal.missionbio.com/datasets/4-cell-lines-AML-multiomics
# Subset of variants was extracted and pathes were extracted and written to 
# demo.h5
# "assays/dna_read_counts/metadata/"
# 'assays/dna_variants/ca/id'
# "assays/dna_read_counts/ra/barcode"
# "assays/dna_variants/ra/barcode"
# "assays/dna_variants/layers/DP"
# "assays/dna_variants/layers/GQ"
# "assays/dna_variants/layers/NGT"
# "assays/dna_variants/layers/AF"
# 'assays/dna_read_counts/ca/id' (Note: This is repeated from the second path.)
# "assays/dna_read_counts/layers/read_counts"
# "/assays/dna_read_counts/ca/CHROM"
# "/assays/dna_read_counts/ca/start_pos"
# "/assays/dna_read_counts/ca/end_pos"


# clusterplot.rds --------------------------------------------------------------
# clusterplot of demo.h5
# h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
# h5 <- h5ToSce(h5_file_path)
# sce <- h5$sce_amp
# se.var <- h5$se_var
# sce <- normalizeReadCounts(sce = sce)
# filteres <- filterVariants(depth.threshold = 10,
#                            genotype.quality.threshold = 30,
#                            vaf.ref = 5, 
#                            vaf.het =  35, 
#                            vaf.hom = 95, 
#                            min.cell = 50,
#                            min.mut.cell = 1,
#                            se.var = se.var,
#                            sce = sce,
#                            shiny = FALSE)
# se.f <- SummarizedExperiment(assays = list(VAF = t(filteres$vaf.matrix.filtered), 
#                                            Genotype = t(filteres$genotype.matrix.filtered),
#                                            Genoqual = t(filteres$genoqual.matrix.filtered)),
#                              rowData = filteres$variant.ids.filtered,
#                              colData = filteres$cells.keep)
# indices_to_keep <- match(filteres$cells.keep, SummarizedExperiment::colData(sce)[[1]], nomatch = 0)
# sce_filtered <- sce[, indices_to_keep]
# SingleCellExperiment::altExp(sce_filtered, "variants") <- se.f
# sce_filtered <- annotateVariants(sce = sce_filtered)
# variants.of.interest <- c("KIT:chr4:55599436:T/C", 
#                           "TET2:chr4:106158216:G/A", 
#                           "FLT3:chr13:28610183:A/G",
#                           "TP53:chr17:7577427:G/A")
# clusterplot <- clusterVariantSelection(sce = sce_filtered,
#                                variants.of.interest = variants.of.interest,
#                                n.clust = 3)

# sce_filtered_demo.rds --------------------------------------------------------
# h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
# h5 <- h5ToSce(h5_file_path)
# sce <- h5$sce_amp
# se.var <- h5$se_var
# sce <- normalizeReadCounts(sce = sce)
# filteres <- filterVariants(depth.threshold = 10,
#                            genotype.quality.threshold = 30,
#                            vaf.ref = 5, 
#                            vaf.het =  35, 
#                            vaf.hom = 95, 
#                            min.cell = 50,
#                            min.mut.cell = 1,
#                            se.var = se.var,
#                            sce = sce,
#                            shiny = FALSE)
# se.f <- SummarizedExperiment(assays = list(VAF = t(filteres$vaf.matrix.filtered), 
#                                            Genotype = t(filteres$genotype.matrix.filtered),
#                                            Genoqual = t(filteres$genoqual.matrix.filtered)),
#                              rowData = filteres$variant.ids.filtered,
#                              colData = filteres$cells.keep)
# indices_to_keep <- match(filteres$cells.keep, SummarizedExperiment::colData(sce)[[1]], nomatch = 0)
# sce_filtered <- sce[, indices_to_keep]
# SingleCellExperiment::altExp(sce_filtered, "variants") <- se.f
# sce_filtered <- annotateVariants(sce = sce_filtered)
# sce_filtered %>% saveRDS('./inst/extdata/sce_filtered_demo.rds')

# sce_norm_demo.rds ------------------------------------------------------------
# Load .h5 file
# h5_file_path <- system.file("extdata", "demo.h5", package = "scafari")
# h5 <- h5ToSce(h5_file_path)
# sce <- h5$sce_amp
# se.var <- h5$se_var
# sce <- normalizeReadCounts(sce = sce)
# saveRDS(sce, './inst/extdata/sce_norm_demo.rds')

# UCSC_hg19_knownCanonical_mock.txt --------------------------------------------
# knownCanonical were downloaded here: 
# https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/
# -> mock was generated with the first 5000 lines (head -n 5000)