#' Function: annotateAmplicons
#' This function takes a SingleCellExperiment object as input and annotates the
#'  stored amplicons.
#'
#' @param sce SingleCellExperiment object containing the single-cell data.
#' @param known.canon Path to jnown canonicals (see vignette)
#' @param shiny If TRUE messages are shown
#'
#' @return A dataframe containing annotated amplicons.
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with a 'counts' assay
#' sce_filtered <- readRDS(system.file("extdata", "sce_filtered_demo.rds",
#'     package = "scafari"
#' ))
#' annotated <- annotateAmplicons(
#'     sce_filtered,
#'     system.file("extdata", "UCSC_hg19_knownCanonical_mock.txt",
#'         package = "scafari"
#'     )
#' )
#' @export
annotateAmplicons <- function(sce, known.canon, shiny = FALSE) {
    # Check that the input is a SingleCellExperiment object
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("`sce` must be a SingleCellExperiment object.")
    }

    # Check if the SCE object has metadata
    if (is.null(sce@metadata)) {
        stop("The SingleCellExperiment object does not contain any metadata.")
    }

    metadata <- sce@metadata

    # Verify that metadata contains the expected element 'genome_version'
    if (is.null(metadata[["genome_version"]])) {
        stop("The metadata must contain 'genome_version'.")
    }

    # Attempt to extract genome version
    genome_version <- tryCatch(
        {
            metadata(sce)[["genome_version"]]
        },
        error = function(e) {
            stop("Failed to extract genome version from metadata: ", e$message)
        }
    )

    sample_name <- tryCatch(
        {
            metadata(sce)[["sample_name"]]
        },
        error = function(e) {
            stop("Failed to extract sample name from metadata: ", e$message)
        }
    )

    if (genome_version == "hg19") {
        # Prepare exon database ------------------------------------------------
        # Read Biomart Exon information and format them
        amps <- as.data.frame(rowData(sce))
        if (shiny) showNotification("Trying to access ensembl...")
        try(mart <- useEnsembl(
            biomart = "ensembl",
            dataset = "hsapiens_gene_ensembl", GRCh = 37
        ))
        if (!exists("mart")) {
            if (shiny) showNotification("Accessing biomart failed.")
            message("Ensemble is not accessible at the moment.")
            amps %>%
                as.data.frame() %>%
                dplyr::mutate(Gene = str_split_i(id, "_", 3)) %>%
                `colnames<-`(c(
                    "Chromosome", "Start", "End", "Amplicon",
                    "Gene"
                )) %>%
                datatable(.,
                    rownames = FALSE, extensions = "Buttons",
                    options = list(
                        pageLength = 10, width = "100%",
                        dom = "Bfrtip",
                        buttons = list(
                            list(extend = "csv", filename = paste0(
                                "scafari_panel_", sample_name
                            )),
                            list(extend = "excel", filename = paste0(
                                "scafari_panel_", sample_name
                            )),
                            list(extend = "pdf", filename = paste0(
                                "scafari_panel_", sample_name
                            )),
                            list(extend = "copy", filename = paste0(
                                "scafari_panel_", sample_name
                            ))
                        )
                    )
                )
        } else {
            get_exon_data <- function(i) {
                exons <- getBM(
                    attributes = c(
                        "ensembl_exon_id",
                        "ensembl_transcript_id_version",
                        "chromosome_name", "exon_chrom_start",
                        "exon_chrom_end", "rank"
                    ),
                    filters = c("chromosome_name", "start", "end"),
                    values = list(
                        gsub("chr", "", amps$seqnames[i]), amps$start[i],
                        amps$end[i]
                    ),
                    mart = mart
                )
                exons$region_id <- amps$id[i] # Include the region ID
                return(exons)
            }

            # Use mclapply for parallel processing
            if (!shiny) {
                all_exon_data <- lapply(seq_len(nrow(amps)), get_exon_data)
            } else {
                withProgress(message = "Processing", value = 0, {
                    total <- nrow(amps)
                    all_exon_data <- lapply(seqlen(total), function(i) {
                        incProgress(1 / total, detail = paste(
                            "Processing amplicon", i,
                            "of", total
                        ))
                        get_exon_data(i)
                    })
                })
            }

            # Combine all exon data into a single data frame
            exon_data <- do.call(rbind, all_exon_data)
            head(exon_data)

            colnames(exon_data) <- c(
                "exon_id", "transcript", "seqnames",
                "start", "end", "rank", "id"
            )
            exon_data$seqnames <- paste0("chr", exon_data$seqnames)

            # Filter out data with invalid annotation values
            exon_data_clean <- exon_data[!startsWith(
                exon_data$exon_id,
                "Error"
            ), ]
            exons.gr <- makeGRangesFromDataFrame(exon_data_clean,
                keep.extra.columns = TRUE
            )

            # Extract canonical transcripts
            canon.path <- known.canon
            known.canon <- read.delim(canon.path,
                header = FALSE,
                col.names = c(
                    "seqnames", "start", "end",
                    "x", "transcript", "ENSG"
                )
            )
            known.canon$transcript <- gsub("\\..*", "", known.canon$transcript)
            exons.gr$transcript <- gsub("\\..*", "", exons.gr$transcript)
            exons.gr.clean <- exons.gr[exons.gr$transcript %in%
                known.canon$transcript, ]
            gene.anno.gr <- makeGRangesFromDataFrame(amps)

            ov <- findOverlaps(gene.anno.gr, exons.gr.clean)

            mcols(gene.anno.gr)["transcript"] <- "-"
            mcols(gene.anno.gr)["Exon"] <- "-"

            gene.anno.gr[queryHits(ov)]$Exon <-
                exons.gr.clean[subjectHits(ov)]$rank
            gene.anno.gr[queryHits(ov)]$transcript <-
                exons.gr.clean[subjectHits(ov)]$transcript

            df <- gene.anno.gr %>%
                as.data.frame() %>%
                tibble::rownames_to_column("id") %>%
                dplyr::mutate(Gene = str_split_i(id, "_", 3)) %>%
                dplyr::select(
                    seqnames, start, end, width, Gene,
                    Exon, transcript
                ) %>%
                `colnames<-`(c(
                    "Chromosome", "Start", "End", "Amplicon length (bp)",
                    "Gene", "Exon", "Canonical Transcript ID"
                )) %>%
                datatable(.,
                    rownames = FALSE, extensions = "Buttons",
                    options = list(
                        pageLength = 10, width = "100%",
                        dom = "Bfrtip",
                        buttons = list(
                            list(
                                extend = "csv",
                                filename = paste0(
                                    "scafari_panel_",
                                    sample_name
                                )
                            ),
                            list(
                                extend = "excel",
                                filename = paste0(
                                    "scafari_panel_",
                                    sample_name
                                )
                            ),
                            list(
                                extend = "pdf",
                                filename = paste0(
                                    "scafari_panel_",
                                    sample_name
                                )
                            ),
                            list(
                                extend = "copy",
                                filename = paste0(
                                    "scafari_panel_",
                                    sample_name
                                )
                            )
                        )
                    )
                )
            return(df)
        }
    } else if (genome_version == "hg38") {
        # MANE annotation
        message("hg38")
        message("MANE annotation is starting. This may take a while.\n")
        mane <- txdbmaker::makeTxDbFromGFF(
            paste0(
                "https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/",
                "release_1.0/MANE.GRCh38.v1.0.ensembl_genomic.gff.gz"
            )
        )
        mane.df <- exonsBy(mane, by = "tx", use.names = TRUE) %>%
            as.data.frame()
        message("Processing MANE\n")
        message("Annotating\n")
        mane.gr <- makeGRangesFromDataFrame(mane.df,
            keep.extra.columns = TRUE,
            na.rm = TRUE
        )

        # Define transcript columns
        gene.anno.gr <- makeGRangesFromDataFrame(rowData(sce))
        mcols(gene.anno.gr)[["Exon"]] <- "-"
        mcols(gene.anno.gr)[["Transcript ID"]] <- "-"
        mcols(gene.anno.gr)[["Gene"]] <- "-"
        ov <- findOverlaps(gene.anno.gr, mane.gr)

        # Until here done
        mcols(gene.anno.gr)[queryHits(ov), ][["Exon"]] <-
            mcols(mane.gr)[subjectHits(ov), ][["exon_rank"]]
        mcols(gene.anno.gr)[queryHits(ov), ]["Transcript ID"] <-
            mcols(mane.gr)[subjectHits(ov), ]["group_name"]
        mcols(gene.anno.gr)["Gene"] <- str_match(
            names(gene.anno.gr),
            "^(.*_)(v\\d_)(.*)_\\d+"
        )[, 4]
        df <- gene.anno.gr %>% as.data.frame()
        df$id <- rownames(df)
        df %>%
            dplyr::select(
                id, seqnames, start, end, width, Gene, Exon,
                `Transcript.ID`
            ) %>%
            `colnames<-`(c(
                "Amplicon ID", "Chromosome", "Start", "End",
                "Amplicon length (bp)", "Gene", "Exon",
                "Canonical Transcript ID"
            )) %>%
            datatable(.,
                extensions = "Buttons",
                options = list(
                    pageLength = 10, width = "100%",
                    dom = "Bfrtip",
                    buttons = list(
                        list(
                            extend = "csv", filename =
                                paste0("scafari_panel_", sample_name)
                        ),
                        list(
                            extend = "excel", filename =
                                paste0("scafari_panel_", sample_name)
                        ),
                        list(
                            extend = "pdf", filename =
                                paste0("scafari_panel_", sample_name)
                        ),
                        list(
                            extend = "copy", filename =
                                paste0("scafari_panel_", sample_name)
                        )
                    )
                ), rownames = FALSE
            )
        return(df)
    } else {
        stop("No proper genome version")
    }
}
