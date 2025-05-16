#' Function: annotateVariants
#' -------------------------------
#' This function takes a SingleCellExperiment object as input and performs variant annotation.
#'
#' @param sce SingleCellExperiment object containing the single-cell data to be annotated.
#' @param shiny A logical flag indicating whether the function is being run in a Shiny application
#' context. Default is FALSE.
#' @param max.var Maximum number of variants to annotate. By default this is 50 to avoid long runtime.
#'
#' @return The function returns an annotated SingleCellExperiment object.
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp()
#' sce_filtered <- readRDS(system.file("extdata", "sce_filtered.rds", package = "scafari"))
#' sce <- annotateVariants(sce_filtered, shiny = FALSE)
#' @export
#'
#' @references https://missionbio.github.io/mosaic/, https://github.com/rachelgriffard/optima
annotateVariants <- function(sce, shiny = FALSE, max.var = 50) {
  # Check that the input is a SingleCellExperiment object
  if (!inherits(sce, "SingleCellExperiment")) {
    stop("`sce` must be a SingleCellExperiment object.")
  }

  if (!inherits(max.var, "numeric")) {
    stop("`max.var` must be a numeric.")
  }

  if (!inherits(shiny, "logical")) {
    stop("`shiny` must be a numeric.")
  }

  # Check if the SCE object has metadata
  if (is.null(metadata(sce))) {
    stop("The SingleCellExperiment object does not contain any metadata.")
  }

  variant.ids.filtered <- rowData(altExp(sce))

  # Check that the input is a SingleCellExperiment object
  if (nrow(altExp(sce)) >= max.var) {
    stop("You try to annotated >= 50 variants. This exceeds `max.var`. If you want to annoate more than 50 variants you need to increase the `max.var` param.")
  }

  metadata <- metadata(sce)

  # Check that genome version is present in metadata
  if (!"genome_version" %in% names(metadata)) {
    stop("The metadata must contain 'genome_version'.")
  }

  if (shiny) {
    message("shiny")
    withProgress(message = "Annotate Variants", value = 0, {
      if (check_MBAPI() == "MissionBio" && metadata[["genome_version"]] == "hg19") {
        if ("annotated" %in% names(metadata(altExp(sce)))) {
          var.vec <- rowData(altExp(sce))$id
        } else {
          var.vec <- rowData(altExp(sce))[["X"]]
        }

        # Bring variants to MissionBio APIs variant format
        matches <- str_match(var.vec, "(chr[0-9XY]+):(\\d+):([ATGC]+)/([ATGC]+)")
        var.mb <- apply(matches, 1, function(x) paste(x[2], x[3], x[4], x[5], sep = "-"))
        variant.ids.filtered.df.anno <- data.frame(
          Gene = character(), Protein = character(), `Coding impact` = character(),
          Function = character(), DANN = character(), ClinVar = character(), dbsnp = character()
        )

        for (var in seq_along(var.mb)) {
          incProgress(1 / length(var.mb), detail = paste0("Annotate variant ", var))

          # Reach MissioBio API
          url <- paste0("https://api.missionbio.io/annotations/v1/variants?ids=", var.mb[var])

          # Parse URL to get variant annotation
          res <- httr::GET(url)
          data <- jsonlite::fromJSON(rawToChar(res$content)) # TODO vary
          annot <- data$annotations %>%
            unlist() %>%
            t() %>%
            as.data.frame()

          if (sum(startsWith(colnames(annot), "function.value")) == 0) {
            annot$function.value <- "-"
          }
          annot <- annot %>%
            mutate("function.value" = tidyr::unite(
              dplyr::select(., starts_with("function.value")),
              col = "function.value",
              sep = "-",
              remove = FALSE
            ) %>%
              pull("function.value")) %>%
            dplyr::select(any_of(c(
              "gene.value", "protein.value", "protein_coding_impact.value",
              "function.value", "impact.value", "clinvar.value", "allele_freq.value", "dbsnp.value"
            ))) %>%
            dplyr::rename(any_of(c(
              Gene = "gene.value",
              Protein = "protein.value",
              `Coding.impact` = "protein_coding_impact.value",
              Function = "function.value",
              DANN = "impact.value",
              ClinVar = "clinvar.value",
              dbsnp = "dbsnp.value",
              `Allele Freq (gnomAD)` = "allele_freq.value"
            ))) %>%
            as.data.frame()
          # rownames(annot) <-  var.mb[var]
          variant.ids.filtered.df.anno <- dplyr::bind_rows(variant.ids.filtered.df.anno, annot) # readRDS('.//input/variant_ids_filtered_df_anno.rds')#
        }


        # Original variant ids as rownames
        rownames(variant.ids.filtered.df.anno) <- variant.ids.filtered$id
        variant.ids.filtered.df.anno$id <- variant.ids.filtered[[1]]
      } else if (check_MBAPI() != "MissionBio" && metadata[["genome_version"]] == "hg19") {
        stop("MissionBio API is not available.")
      } else if (metadata[["genome_version"]] == "hg38") {
        snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

        if ("annotated" %in% names(metadata(altExp(sce)))) {
          variant.ids.filtered <- rowData(altExp(sce))$id
        } else {
          variant.ids.filtered <- rowData(altExp(sce))[["X"]]
        }


        anno <- data.frame(
          "refsnp_source" = c(), "refsnp_id" = c(),
          "chr_name" = c(), "chrom_start" = c(),
          "chrom_end" = c(), "consequence_type_tv" = c(),
          "clinical_significance" = c(), "ENSEMBL" = c(), "id" = c()
        )

        for (var in variant.ids.filtered[1:3]) {
          message(paste0("annotate ", var))
          var.tmp <- str_match(var, "([0-9XY]+):(\\d+)") %>%
            {
              paste(.[2], .[3], .[3], sep = ":")
            }
          tryCatch(
            {
              anno.tmp <- getBM(
                attributes = c(
                  "refsnp_source", "refsnp_id", "chr_name", "chrom_start", 
                  "chrom_end",
                  "consequence_type_tv", "clinical_significance", "ENSEMBL"
                ),
                filters = "chromosomal_region",
                values = var.tmp,
                mart = snpmart
              )
              anno.tmp$ID <- var
            },
            error = function(e) {
              # Return NA for this SNP if an error occurs
              data.frame(
                refsnp_source = NA, refsnp_ID = NA, chr_name = NA, 
                chrom_start = NA,
                chrom_end = NA, consequence_type_tv = NA, 
                clinical_significance = NA,
                ensembl_gene_name = NA
              )
            }
          )
          anno <- rbind(anno, anno.tmp)
        }

        lookup <- select(org.Hs.eg.db, keys = unique(anno$ensembl_gene_name), 
                         keytype = "ENSEMBL", columns = "SYMBOL") 
        variant.ids.filtered.df.anno <- merge(anno, lookup, all = TRUE) %>%
          dplyr::select(any_of(c(
            "ID", "SYMBOL", "Position", "consequence_type_tv",
            "clinical_significance",
            "refsnp_id"
          ))) %>%
          dplyr::rename(any_of(c(
            Gene = "SYMBOL",
            `chr_name` = "Chr",
            `Consequence` = "consequence_type_tv",
            `Clinical significance` = "clinical_significance",
            dbsnp = "refsnp_id"
          ))) %>%
          tidyr::separate(ID, into = c("Chr", "Start", "Ref", "Alt"), 
                          sep = ":", remove = FALSE) %>%
          as.data.frame()
      } else {
        message("Issue with MissionBio API. Using Biomart")
      }
    })


    # Without shiny logic
  } else {
    if (check_MBAPI() == "MissionBio" && 
        metadata[["genome_version"]] == "hg19") {
      if ("annotated" %in% names(metadata(altExp(sce)))) {
        var.vec <- rowData(altExp(sce))$id
      } else {
        var.vec <- rowData(altExp(sce))[["X"]]
      }

      # Bring variants to MissionBio APIs variant format
      matches <- str_match(var.vec, "(chr[0-9XY]+):(\\d+):([ATGC]+)/([ATGC]+)")
      var.mb <- apply(matches, 1, function(x) paste(x[2], x[3], x[4], x[5], 
                                                    sep = "-"))
      variant.ids.filtered.df.anno <- data.frame(
        Gene = character(), Protein = character(), 
        `Coding impact` = character(),
        Function = character(), DANN = character(), 
        ClinVar = character(), dbsnp = character()
      )

      for (var in seq_along(var.mb)) {
        # Reach MissioBio API
        url <- paste0("https://api.missionbio.io/annotations/v1/variants?ids=", 
                      var.mb[var])

        # Parse URL to get variant annotation
        res <- httr::GET(url)
        data <- jsonlite::fromJSON(rawToChar(res$content)) # TODO vary
        annot <- data$annotations %>%
          unlist() %>%
          t() %>%
          as.data.frame()

        if (sum(startsWith(colnames(annot), "function.value")) == 0) {
          annot$function.value <- "-"
        }
        annot <- annot %>%
          mutate("function.value" = tidyr::unite(
            dplyr::select(., starts_with("function.value")),
            col = "function.value",
            sep = "-",
            remove = FALSE
          ) %>%
            pull("function.value")) %>%
          dplyr::select(any_of(c(
            "gene.value", "protein.value", "protein_coding_impact.value",
            "function.value", "impact.value", "clinvar.value", 
            "allele_freq.value", "dbsnp.value"
          ))) %>%
          dplyr::rename(any_of(c(
            Gene = "gene.value",
            Protein = "protein.value",
            `Coding.impact` = "protein_coding_impact.value",
            Function = "function.value",
            DANN = "impact.value",
            ClinVar = "clinvar.value",
            dbsnp = "dbsnp.value",
            `Allele Freq (gnomAD)` = "allele_freq.value"
          ))) %>%
          as.data.frame()
        # rownames(annot) <-  var.mb[var]
        variant.ids.filtered.df.anno <- dplyr::bind_rows(variant.ids.filtered.df.anno, annot) # readRDS('.//input/variant_ids_filtered_df_anno.rds')#
      }


      # Original variant ids as rownames
      rownames(variant.ids.filtered.df.anno) <- variant.ids.filtered$id
      variant.ids.filtered.df.anno$id <- variant.ids.filtered[[1]]
    } else if (check_MBAPI() != "MissionBio" && metadata[["genome_version"]] == "hg19") {
      # TODO
    } else if (metadata[["genome_version"]] == "hg38") {
      snpmart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

      anno <- data.frame(
        "refsnp_source" = c(), "refsnp_id" = c(),
        "chr_name" = c(), "chrom_start" = c(),
        "chrom_end" = c(), "consequence_type_tv" = c(),
        "clinical_significance" = c(), "ensembl_gene_name" = c(), "id" = c()
      )

      for (var in variant.ids.filtered[1:3]) {
        message(paste0("annotate ", var))
        var.tmp <- str_match(var, "([0-9XY]+):(\\d+)") %>%
          {
            paste(.[2], .[3], .[3], sep = ":")
          }
        tryCatch(
          {
            anno.tmp <- getBM(
              attributes = c(
                "refsnp_source", "refsnp_id", "chr_name", "chrom_start", 
                "chrom_end","consequence_type_tv", "clinical_significance", 
                "ensembl_gene_name"
              ),
              filters = "chromosomal_region",
              values = var.tmp,
              mart = snpmart
            )
            anno.tmp$ID <- var
          },
          error = function(e) {
            # Return NA for this SNP if an error occurs
            data.frame(
              refsnp_source = NA, refsnp_ID = NA, chr_name = NA, 
              chrom_start = NA, chrom_end = NA, consequence_type_tv = NA, 
              clinical_significance = NA, ensembl_gene_name = NA
            )
          }
        )
        anno <- rbind(anno, anno.tmp)
      }

      lookup <- select(org.Hs.eg.db, keys = unique(anno$ensembl_gene_name), 
                       keytype = "ENSEMBL", columns = "SYMBOL") 
      variant.ids.filtered.df.anno <- merge(anno, lookup, all = TRUE) %>%
        dplyr::select(any_of(c(
          "ID", "SYMBOL", "Position", "consequence_type_tv",
          "clinical_significance",
          "refsnp_id"
        ))) %>%
        dplyr::rename(any_of(c(
          Gene = "SYMBOL",
          `chr_name` = "Chr",
          `Consequence` = "consequence_type_tv",
          `Clinical significance` = "clinical_significance",
          dbsnp = "refsnp_id"
        ))) %>%
        tidyr::separate(ID, into = c("Chr", "Start", "Ref", "Alt"), 
                        sep = ":", remove = FALSE) %>%
        as.data.frame()
    } else {
      message("Issue with MissionBio API. Using Biomart")
    }
  }
  rowData(altExp(sce)) <- variant.ids.filtered.df.anno
  metadata(altExp(sce))[["annotated"]] <- TRUE
  return(sce)
}
