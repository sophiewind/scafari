# Check if MissionBio is accessible
check_MBAPI <- function() {
  tryCatch({
    response <- httr::GET('https://api.missionbio.io/annotations/v1/variants?ids=17:7578211:7578211')
    status <- httr::status_code(response)
    if (status >= 200 && status < 300) {
      data = jsonlite::fromJSON(rawToChar(response$content))
      if (length(data) > 2){
        return('MissionBio')
      }
    } else {
      return('Biomart')
      stop(paste("URL not available. Status code:", status))
    }
  }, error = function(e) {
    return('Biomart')
    stop(paste("Error accessing URL:", e$message))
  })
}

## Check h5
check_h5_structure <- function(h5_file) {
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
  full_paths <- paste(h5_structure$group, h5_structure$name, sep="/")
  missing_paths <- character(0)
  for (path in paths_to_check) {
    shinyjs::html("text", paste0('Checking ', path, '<br>'), add = F)
    
    
    if (!(path %in% full_paths)) {
      missing_paths <- c(missing_paths, path)
      print()
      shinyjs::html("text",paste0(path, ' is missing!<br>'), add = F)
      
    }
    Sys.sleep(0.01)
  }
  if (length(missing_paths) > 0) {
    shinyjs::html("text", paste0("<p style='color:red'>", missing_paths, " missing."))
    return(list(valid = FALSE, missing = missing_paths))
  } else {
    shinyjs::html("text", paste0("<p style='color:green'>", "All pathes are available."))
  }
  return(list(valid = TRUE))
}
