# STAT -------------------------------------------------------------------------
#' Default shiny theme
#'
#' @return gg color scheme
#'
#' @keywords internal
theme_shiny <- function() {
    theme_minimal() +
        theme(
            panel.grid = element_blank(),
            title = element_text(size = 20),
            text = element_text(size = 16)
        )
}

#' Default ggplot scheme
#'
#' @return gg color scheme
#'
#' @keywords internal
theme_default <- function() {
    theme_minimal() +
        theme(
            panel.grid = element_blank(),
            title = element_text(size = 12),
            text = element_text(size = 10)
        )
}

#' Default ggplot colors
#'
#' @return gg colors for plotting
#'
#' @keywords internal
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq(1, n)]
}

mycols <- c(
    "WT" = "#414487FF", "Hom" = "#D44292", "Het" = "#F6A97A",
    "Missing" = "#868686FF"
)
mycols.ngt <- c(
    `WT` = "#414487FF", `Hom` = "#D44292", `Het` = "#F6A97A",
    `Missing` = "#868686FF"
)



colors_vaf <- circlize::colorRamp2(
    c(0, 50, 100),
    c("#414487FF", "#F6A97A", "#D44292")
)


# Chromosome color palette
chromosomes <- c(paste0("chr", seq(1, 21)), "chrX", "chrY")

# Get colors from the viridis "magma" palette
colors <- gg_color_hue(length(chromosomes))

# Create a named list for the color palette
chr_palette <- setNames(colors, chromosomes)

# Upload -----------------------------------------------------------------------
#' Check h5
#' @param h5_file and path to the h5 file
#'
#' @return boolean if h5 is valid
#'
#' @keywords internal
checkH5 <- function(h5_file) {
    paths_to_check <- c(
        "/assays/dna_variants/layers/DP",
        "/assays/dna_variants/layers/GQ",
        "/assays/dna_variants/layers/NGT",
        "/assays/dna_variants/layers/AF",
        "/assays/dna_read_counts/ca/id",
        "/assays/dna_read_counts/layers/read_counts",
        "/assays/dna_read_counts/ca/CHROM",
        "/assays/dna_read_counts/ca/start_pos",
        "/assays/dna_read_counts/ca/end_pos"
    )
    h5_structure <- h5ls(h5_file)
    full_paths <- paste(h5_structure$group, h5_structure$name, sep = "/")
    missing_paths <- character(0)
    for (path in paths_to_check) {
        shinyjs::html("text", paste0("Checking ", path, "<br>"), add = FALSE)


        if (!(path %in% full_paths)) {
            missing_paths <- c(missing_paths, path)
            shinyjs::html("text", paste0(path, " is missing!<br>"), add = FALSE)
        }
        Sys.sleep(0.01)
    }
    if (length(missing_paths) > 0) {
        shinyjs::html("text", paste0(
            "<p style='color:red'>", missing_paths,
            " missing."
        ))
        return(list(valid = FALSE, missing = missing_paths))
    } else {
        shinyjs::html("text", paste0(
            "<p style='color:green'>",
            "All pathes are available."
        ))
    }
    return(list(valid = TRUE))
}



#' Check if MissionBio is accessible
#'
#' @return boolean if MissionBio API is currently available
#'
#' @keywords internal
check_MBAPI <- function() {
    response <-
        httr::GET(paste0(
            "https://api.missionbio.io/annotations/v1/",
            "variants?ids=17:7578211:7578211"
        ))
    status <- httr::status_code(response)
    if (status >= 200 && status < 300) {
        data <- jsonlite::fromJSON(rawToChar(response$content))
        if (length(data) > 2) {
            return("MissionBio")
        } else {
            return("Biomart")
        }
    } else {
        stop("URL not available. Status code:", status)
    }
}


#' CheckSce 
#' 
#' Checks sce objects.
#' 
#' @return TRUE if input is valid; otherwise an error is thrown
#' 
#' @keywords internal
#' @noRd
checkSce <- function(sce, variants = FALSE, gt = FALSE){
    # Check that the input is a SingleCellExperiment object
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("The input must be a SingleCellExperiment object.")
    }
    
    if (variants){
        # Check for the presence of 'variants' altExp and required assays
        if (!"variants" %in% altExpNames(sce)) {
            stop("The SingleCellExperiment object must contain 'variants' ",
                "as an alternate experiment.")
        }
        if (!"VAF" %in% assayNames(altExp(sce, "variants"))) {
            stop("The 'variants' alternate experiment must contain a 'VAF'",
            "assay.")
        }
    }
    
    if (gt){
        if (!all(c("VAF", "Genotype") %in% 
                assayNames(altExp(sce, "variants")))) {
            stop("The 'variants' alternate experiment must contain 'VAF' ",
                "and Genotype' assays.")
        }
    }
    
    TRUE
}

