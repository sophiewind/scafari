# Load necessary packages
library(shiny)
library(shinycssloaders)
library(DT)
library(dplyr)
library(waiter)
library(ggplot2)
library(GenomicRanges)
library(tibble)
library(rhdf5)
library(ComplexHeatmap)
library(karyoploteR)
library(stringr)
library(reshape2)
library(biomaRt)  # NEW
library(clusterProfiler)  # NEW
library(org.Hs.eg.db)  # NEW
library(shinyjs) # NEW
library(shinyBS) # NEW
library(shinycustomloader) # NEW and un
library(factoextra)
library(markdown)
library(tidyr)

h5closeAll()  

# Define the UI
ui <- navbarPage(
  title = 'scafari',
  
  # Upload tab -----------------------------------------------------------------
  tabPanel("Upload",
           fluidPage(
             use_waiter(),  # Include waiter in the UI
             useShinyjs(),  # Enable shinyjs
             width = 3,  # Set sidebar width to 3 columns (out of 12)
             h3("Data Upload"),
             fileInput("upload", "Upload HDF5 File", accept = '.h5'),
             verbatimTextOutput("class_output"),
             verbatimTextOutput("file_contents"),  # Output area for file contents  
             div(id = "text"),
             tableOutput("files"))),
  
  # Sequencing tab -------------------------------------------------------------
  tabPanel("Sequencing",
           fluidPage(
             use_waiter(),  # Include waiter in the UI
             useShinyjs(),  # Enable shinyjs
             # Conditional display based on file upload
             conditionalPanel(
               condition = "output.file_ready",  # Check if a file is uploaded
               fluidRow(
                 h2('Overview'),
                 column(width = 6,
                        plotOutput("seq_plot1", height = "300px"),
                        withLoader(plotOutput("seq_plot2", height = "300px"), loader = 'dnaspin')
                 ),
                 column(width = 6,
                        withLoader(plotOutput("seq_plot3", height = "600px"), loader = 'dnaspin')
                 )
               ),
               hr(),
               fluidRow(
                 h2('Sample information'),
                 withLoader(dataTableOutput("data_table_sample"), loader = 'dnaspin')
               ),
               hr(),
               fluidRow(
                 h2('Sequencing information'),
                 withLoader(dataTableOutput("data_table_sequencing"), loader = 'dnaspin')
               ),
               hr(),
               fluidRow(
                 h2('Mapping'),
                 withLoader(dataTableOutput("data_table_mapping"), loader = 'dnaspin')
               )
             ),
             # Message displayed when the file is not uploaded
             conditionalPanel(
               condition = "!output.file_ready",  # Check if no file is uploaded
               fluidRow(
                 style = "text-align: center; margin-top: 50px;", 
                 p(style = "color: black;", "Please upload a file to view sequencing information.")                   )
             ),
             hr(),
           )
  ),
  
  # Panel tab ------------------------------------------------------------------
  tabPanel("Panel",
           fluidPage(
             use_waiter(),  # Include waiter in the UI
             useShinyjs(),  # Enable shinyjs
             conditionalPanel(
               condition = "output.file_ready",
               h2("Panel Information"),  # Header for the panel
               fluidRow(
                 column(width=6,
                        withLoader(plotOutput("panel_plot1")), loader = 'dnaspin'),
                 column(width=6,
                        withLoader(plotOutput("panel_plot2")), loader = 'dnaspin')
               ),
               fluidRow(
                 h2('Amplicon Overview'),
                 withLoader(dataTableOutput("data_table_overview"))
               ),
               fluidRow(
                 h2('Amplicon Distribution'),
                 withLoader(plotOutput("panel_plot3", height = "800px"), loader = 'dnaspin')
               ),
               fluidRow(
                 h2('Amplicon Performance'),
                 withLoader(plotOutput("panel_plot4"), loader = 'dnaspin')
               ),
               fluidRow(
                 h2('Panel uniformity plot',bsButton("pu", label = "", icon = icon("info"), style = "info", size = "extra-small")),
                 bsPopover(
                   id = "pu",
                   title = "Panel uniformity",
                   content = paste0(
                     "The plot indicates the normalized read counts per amplicon. Amplicons that meet at least 20% of the average depth of coverage are plotted in blue others in orange."
                   ),
                   placement = "right",
                   trigger = "hover",
                   options = list(container = "body")
                 ),
                 withLoader(plotOutput("panel_plot5"), loader = 'dnaspin')
               )
             ),
             # Message displayed when the file is not uploaded
             conditionalPanel(
               condition = "!output.file_ready",  # Check if no file is uploaded
               fluidRow(
                 style = "text-align: center; margin-top: 50px;", 
                 p(style = "color: black;", "Please upload a file to view sequencing information.")                   )
             ),
             hr(),
          )
  ),
  
  # Variants tab ---------------------------------------------------------------
  tabPanel("Variants",
           fluidPage(
             use_waiter(),  # Include waiter in the UI
             useShinyjs(),  # Enable shinyjs
             wellPanel(fluidRow(
               h2('Filtering Parameters'),
               column(width=6,
                      style = "background-color: #f7f7f7; padding: 10px; margin-bottom: 20px;",
                      
                      numericInput("depth_threshold", "Depth Threshold:", value = 10, min = 1),
                      bsPopover("depth_threshold", title = "Depth Threshold",
                                content = "Enter the minimum depth of coverage for filtering.", 
                                trigger = "hover"),
                      
                      numericInput("genotype_quality_threshold", "Genotype Quality Threshold:", value = 30, min = 1),
                      bsPopover("genotype_quality_threshold", title = "Genotype Quality Threshold",
                                content = "Enter the minimum genotype quality score.", 
                                trigger = "hover"),
                      
                      numericInput("vaf_ref", "VAF Reference Threshold:", value = 5, min = 0, max = 100, step = 10),
                      bsPopover("vaf_ref", title = "VAF Reference Threshold",
                                content = "Define the upper VAF threshold for wildtype variants", 
                                trigger = "hover"),
                      
                      numericInput("vaf_hom", "VAF Homozygous Threshold:", value = 95, min = 0, max = 100, step = 10),
                      bsPopover("vaf_hom", title = "VAF Homozygous Threshold",
                                content = "Define the lower VAF threshold for homozygous variants.", 
                                trigger = "hover")
               ),
               column(width=6,
                      style = "background-color: #f7f7f7; padding: 10px; margin-bottom: 20px;",
                      numericInput("vaf_het", "VAF Heterozygous Threshold:", value = 35, min = 0, max = 100, step = 10),
                      bsPopover("vaf_het", title = "VAF Heterozygous Threshold",
                                content = "Define the lower VAF threshold for heterozygous variants.", 
                                trigger = "hover"),
                      
                      numericInput("min_cell_pt", "Minimum Cell Percentage:", value = 50, min = 0, max = 100),
                      bsPopover("min_cell_pt", title = "Minimum Cell Percentage",
                                content = "Enter the minimum percentage of cells with the variant (%).", 
                                trigger = "hover"),
                      
                      numericInput("min_mut_cell_pt", "Minimum Mutated Cell Percentage:", value = 1, min = 0, max = 100),
                      bsPopover("min_mut_cell_pt", title = "Minimum Mutated Cell Percentage",
                                content = "Enter the minimum percentage of mutated cells required (%).", 
                                trigger = "hover"),
               ),
               actionButton("filter_btn", "Apply Filtering", icon = icon("filter"), 
                            style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
             ),
             conditionalPanel(
               condition = "output.plots_visible == true",
               h2("Variant Information"),
               fluidRow(
                 column(width=6,
                        withLoader(plotOutput("var_plot1", height = "300px")), loader = 'dnaspin'),
                 column(width=6,
                        withLoader(plotOutput("var_plot2", height = "300px")), loader = 'dnaspin')
               ),
               
               fluidRow(
                 h2('Overview Filtered Variants'),
                 dataTableOutput("data_table_var")  
               ),
               fluidRow(
                 h2('Variant Allele Frequency'),
                 column(11,
                        withLoader(plotOutput("var_plot3", height = "600px"), loader = 'dnaspin')
                 ),
                 column(1,
                        plotOutput("legend")),
               ),
               fluidRow(
                 h2('Genotype of Filtered Variants'),
                 column(width=6,
                        withLoader(plotOutput("var_plot4")), loader = 'dnaspin'),
                 column(width=6,
                        withLoader(plotOutput("var_plot5")), loader = 'dnaspin')
               ),
             ),
             hr(),
           )
  ),
  
  # Explore Variants tab -------------------------------------------------------
  tabPanel("Explore profiles",
           fluidPage(
             use_waiter(),  # Include waiter in the UI
             useShinyjs(),  # Enable shinyjs
             fluidRow(
               sidebarLayout(
                 sidebarPanel(width = 3,
                              tags$style(HTML(".dataTable {font-size: 12px; overflow-y: scroll; }")),
                              dataTableOutput("data_table_var2"),  
                              hr(),
                              actionButton('submit_var', label = 'Select variants', icon = icon('play-circle'), 
                                           style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),),
                 
                 mainPanel(h1('Explore profiles'),
                           plotOutput('hm_1', height = "800px" ),
                           
                 )
               ),
               div(style = "display:inline-block; float:center", 
                   actionButton('continue_var', 'Continue with this variant selection', icon = icon('thumbs-up'), 
                                style="color: #fff; background-color: #337ab7; border-color: #2e6da4")),
               hr()
             ),
             conditionalPanel(
               condition = "output.continue == true",
               h1('Identify cell clusters'),
               h2('Select number of clusters'),
               withLoader(plotOutput("kneeplot")), loader = 'dnaspin',
               fluidRow(
                 numericInput("n_clust", "Number of clusters:", value = 3, min = 2, max = 10, step = 1), 
                 
                 div(style = "display:inline-block; float:center", 
                     actionButton('kmeans_btn', label = 'Start kmeans', icon = icon('play'), 
                                  style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                 ),
                 hr()
               ),
               conditionalPanel(
                 condition = "output.plots_visible_2 == true",
                 
                 h2('Cluster cells by VAF'),
                 h4('Variants included in Clustering:'),
                 uiOutput("selected_rows_2"),
                 fluidRow(withLoader(plotOutput("cluster_plot", height = "800px"), loader = 'dnaspin')),
                 
                 h2('\n'),
                 fluidRow(withLoader(plotOutput("cnv_plot3", height = "800px"), loader = 'dnaspin')),
                 hr(),
                 h2('Explore variant profiles'),
                 h3('VAF in clusters'),
                 fluidRow(withLoader(plotOutput("cnv_plot4"), loader = 'dnaspin')),
                 
                 h3('VAF distribution map'),
                 fluidRow(withLoader(plotOutput("cnv_plot6"), loader = 'dnaspin')),
                                  
                 h3('Numerical genotype in clusters'),
                 fluidRow(plotOutput('ana_bar')),
               )
             ) 
           ),
           hr(),
  )
)

# Server -----------------------------------------------------------------------
server <- function(input, output, session) {

  # Setup statics --------------------------------------------------------------
  # Define plotting statics ---
  gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  
  mycols <- c('WT' = "#414487FF", 'Hom' = "#D44292", 'Het' = "#F6A97A",
              'Missing' = "#868686FF")
  mycols.ngt <- c(`WT` = "#414487FF", `Hom` = "#D44292", `Het` = "#F6A97A",
                  `Missing` = "#868686FF")
  
  theme_default <- function() {
    theme_minimal() +
      theme(
        panel.grid = element_blank(),
        title = element_text(size = 20),
        text = element_text(size = 16)
      )
  }
  
  # Chromosome color palette
  chromosomes <- c(paste0("chr", 1:21), "chrX", "chrY")
  
  # Get colors from the viridis "magma" palette
  colors <- gg_color_hue(length(chromosomes))
  
  # Create a named list for the color palette
  chr_palette <- setNames(colors, chromosomes)
  # ---
  
  # Setup reactivity ---
  plots_visible <- reactiveVal(FALSE)
  plots_visible_2 <- reactiveVal(FALSE)
  continue <- reactiveVal(FALSE)
  
  shinyjs::disable("filter_btn")
  shinyjs::disable("kmeans_btn")
  
  # Upload settings
  options(shiny.maxRequestSize=500*1024^2) 
  # ---
  
  # End of setting up  ---------------------------------------------------------
  source('functions.R')
  
  # Process input --------------------------------------------------------------
  ## Check input ----------------------------------------------------------------
  h5_file <- reactive({
    req(input$upload)  # Ensure input$upload exists
    
    # Open the HDF5 file
    h5_file <- H5Fopen(input$upload$datapath)
    
    # Check existence of important slots in input
    check_result <- check_h5_structure(h5_file)
    
    # If input is valid upload gets deactivated and buttons are activated
    if (check_result$valid) {
      shinyjs::enable("filter_btn")  # Enable button if valid
      shinyjs::disable('upload')  # Deactivate input if valid
      output$status <- renderText("File uploaded successfully!")
      return(h5_file)
    } else {
      showNotification(
        paste("Error: The following paths are missing in the H5 file:", 
              paste(check_result$missing, collapse = ", ")),
        type = "error",
        duration = NULL
      )
      return(NULL)  # Return NULL if invalid
    }
  })
  
  # Display conditional panels if file is uploaded successfully
  output$file_ready <- reactive({
    !is.null(h5_file())
  })
  outputOptions(output, "file_ready", suspendWhenHidden = FALSE) 
  
  
  ## Read input -----------------------------------------------------------------
  # Read metadata
  metadata <- reactive({
    req(h5_file())
    return(as.data.frame(unlist(h5read(h5_file(),"assays/dna_read_counts/metadata/"))))
  })
  
  # Read variant ids
  variant.ids <- reactive({
    req(h5_file())
    h5read(h5_file(), "assays/dna_variants/ca/id")
  })
  
  variant.ids.anno <- reactive({
    NULL
  })
  
  # Read depth matrix
  depth_matrix <- reactive({
    req(h5_file())
    depth_data <- h5read(h5_file(), "assays/dna_variants/layers/DP")
    return(as.data.frame(t(unlist(depth_data))))
  })
  
  # Read genotype quality matrix
  genoqual_matrix <- reactive({
    req(h5_file())
    genoqual_data <- h5read(h5_file(), "assays/dna_variants/layers/GQ")
    t(as.data.frame(unlist(genoqual_data)))
  })
  
  # Read genotype matrix
  genotype_matrix <- reactive({
    req(h5_file())
    genotype_data <- h5read(h5_file(), "assays/dna_variants/layers/NGT")
    t(as.data.frame(unlist(genotype_data)))
  })
  
  # Read variant allele frequency matrix
  vaf_matrix <- reactive({
    req(h5_file())
    vaf_data <- h5read(h5_file(), "assays/dna_variants/layers/AF")
    as.data.frame(t(unlist(vaf_data)))
  })
  
  # Read amplicons
  amplicons <- reactive({
    req(h5_file())
    h5read(h5_file(), "assays/dna_read_counts/ca/id")
  })
  
  # Read read counts
  read.counts.df <- reactive({
    req(h5_file())
    read_counts_data <- h5read(h5_file(), "assays/dna_read_counts/layers/read_counts")
    read_counts_df <- as.data.frame(t(read_counts_data))
    colnames(read_counts_df) <- amplicons()  # Set column names to amplicons
    return(read_counts_df)
  })
  
  # Read gene annotations
  gene.anno.df <- reactive({
    req(h5_file())
    df <- data.frame(
      seqnames = h5read(h5_file(), "/assays/dna_read_counts/ca/CHROM"),
      start = h5read(h5_file(), "/assays/dna_read_counts/ca/start_pos"),
      end = h5read(h5_file(), "/assays/dna_read_counts/ca/end_pos"),
      id = h5read(h5_file(), "/assays/dna_read_counts/ca/id")
    )
    
    # Check if seqnames starts with 'chr' and prepend 'chr' if not
    df$seqnames <- ifelse(startsWith(df$seqnames, 'chr'), 
                          df$seqnames, 
                          paste0('chr', df$seqnames))
    
    return(df)
  })
  
  # Panel analyses--------------------------------------------------------------

  
  ## Normalize read counts -----------------------------------------------------
  # Calculate row sum threshold
  # Set threshold to 1/10 of the 10th smallest total read count
  rowsum.threshold <- reactive({
    req(read.counts.df())
    print(mean(as.matrix(read.counts.df())))
    
    sort(rowSums(read.counts.df()), decreasing = TRUE)[10] / 10 
  })
  
  # Extract cells that will not be removed due to small total read counts
  keep.tf <- reactive({
    req(read.counts.df(), rowsum.threshold())
    rowSums(read.counts.df()) > rowsum.threshold() | rowSums(read.counts.df()) >  40*as.numeric(metadata()['n_amplicons',])
  })
  
  # Normalize read counts
  read.counts.df.norm <- reactive({
    withProgress(message = 'Normalize read counts', value = 0, {
      incProgress(0.1)
      req(read.counts.df(), keep.tf())
      
      # Scale read counts by mean read counts per cell
      read.counts.df.tmp <- read.counts.df() / (apply(read.counts.df(), 1, mean) + 1)  # Scale values, adding 1 to avoid division by zero
      incProgress(0.3)
      
      # Normalize over cells using median reads per amp in cells kept
      read.counts.df.tmp <- t(t(read.counts.df.tmp) / (apply(read.counts.df.tmp[keep.tf(), ], 2, median) + 0.05))  # Normalize by median
      incProgress(0.3)
      
      read.counts.df.tmp <- read.counts.df.tmp * 2  # Scale the normalized counts by a factor of 2
      incProgress(0.3)
      
      as.data.frame(read.counts.df.tmp)  # Convert back to data frame
    })
  })
  
  # Variant analyses -----------------------------------------------------------
  ## Variant filtering ---------------------------------------------------------
  # Filter variants after filter button hit  
  observeEvent(input$filter_btn, {
    plots_visible(TRUE) 
    read.counts.df.norm <- read.counts.df.norm()
    withProgress(message = 'Filter variants', value = 0, {
      incProgress(0, detail = "Depth Filtering...")
      
      # Initial filtering flags
      dp_tf <- depth_matrix() < input$depth_threshold
      incProgress(1/5, detail = "Genotype Quality Filtering...")
      gq_tf <- genoqual_matrix() < input$genotype_quality_threshold
      incProgress(1/5, detail = "Variant Allel Frequency and Numerical Genotype Filtering...")
      vaf_ref_tf <- (vaf_matrix() > input$vaf_ref) & (genotype_matrix() == 0)
      vaf_hom_tf <- (vaf_matrix() < input$vaf_hom) & (genotype_matrix() == 2)
      vaf_het_tf <- (vaf_matrix() < input$vaf_het) & (genotype_matrix() == 1)
      
      # First pass: Determine which entries to keep
      keep <- !(dp_tf | gq_tf | vaf_ref_tf | vaf_hom_tf | vaf_het_tf)
      genotype.matrix <- genotype_matrix()
      vaf.matrix <- vaf_matrix()
      genotype.matrix[!keep] <- 3
      vaf.matrix[genotype.matrix == 3] <- -1
      
      # Filter based on cell and mutation counts
      incProgress(1/5, detail = "Cell and Variant Filtering...")
      num_cells <- nrow(genotype.matrix)
      num_variants <- ncol(genotype.matrix)
      
      cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) {x %in% 0:2})) > num_cells * input$min_cell_pt / 100
      mut_cell_num_keep_tf <- colSums(apply(genotype.matrix, 2, function(x) { x %in% 1:2 })) > num_cells * input$min_mut_cell_pt / 100
      variant_keep_tf <- cell_num_keep_tf & mut_cell_num_keep_tf  
      
      # Second pass filtering
      v_names <- variant.ids()
      filtered_variant_names <- v_names[variant_keep_tf]
      cell_variants_keep_tf <- rowSums(genotype.matrix != 3) > num_variants * input$min_cell_pt / 100
      vaf_matrix_filtered <-vaf.matrix[cell_variants_keep_tf, variant_keep_tf]  
      
      incProgress(1/5, detail = "Processing Filtered Cells and Variants")
      genotype_matrix_filtered <- genotype.matrix[cell_variants_keep_tf, variant_keep_tf]  
      genoqual_matrix_filtered <- genoqual_matrix()[cell_variants_keep_tf, variant_keep_tf]
      read.counts.df.norm <- read.counts.df.norm[cell_variants_keep_tf,]
      
      variant.ids.filtered <- (v_names[variant_keep_tf]) ## Achtung NA
      variant.ids.filtered <- variant.ids.filtered[!is.na(variant.ids.filtered)]
      rownames(vaf_matrix_filtered) <- NULL
    })
    
    ## Annotate variants -------------------------------------------------------
    # Test if missionbio api is available
    withProgress(message = 'Annotate Variants', value = 0, {
      
      if (check_MBAPI() == 'MissionBio' && metadata()['genome_version',] == 'hg19'){
        # Bring variants to MissionBio APIs variant format
        var.mb <- apply(variant.ids.filtered, 1, function(x){str_match(x, "(chr[0-9XY]+):(\\d+):([ATGC]+)/([ATGC]+)") %>%
            { paste(.[2], .[3], .[4], .[5], sep = "-") }})
        variant.ids.filtered.df.anno <- data.frame(Gene = character(), Protein = character(), `Coding impact` = character(),
                                                   Function = character(), DANN = character(), ClinVar = character() , dbsnp = character())
        
        for (var in seq_along(var.mb)) {
          incProgress(1/length(var.mb), detail = paste0('Annotate variant ', var))
          
          # Reach MissioBio API
          url <- paste0("https://api.missionbio.io/annotations/v1/variants?ids=", var.mb[var])
          
          # Parse URL to get variant annotation
          res = httr::GET(url)
          data = jsonlite::fromJSON(rawToChar(res$content))  # TODO vary
          annot <- data$annotations %>%
            unlist() %>%
            t() %>%
            as.data.frame() 
          
          if (sum (startsWith(colnames(annot), 'function.value')) == 0){
            annot$function.value <- '-'
            
          }
          annot <- annot %>%
            mutate('function.value' = tidyr::unite(
            dplyr::select(., starts_with("function.value")),
            col = 'function.value',
            sep = "-",
            remove = FALSE
          ) %>%
            pull('function.value')) %>% 
            dplyr::select(any_of(c('gene.value', 'protein.value', 'protein_coding_impact.value',
                                   'function.value', 'impact.value', 'clinvar.value', "allele_freq.value", 'dbsnp.value'))) %>%
            dplyr::rename(any_of(c(
              Gene = 'gene.value',
              Protein = 'protein.value',
              `Coding.impact` = 'protein_coding_impact.value',
              Function = 'function.value',
              DANN = 'impact.value',
              ClinVar = 'clinvar.value',
              dbsnp = 'dbsnp.value',
              `Allele Freq (gnomAD)` = "allele_freq.value"))) %>%
            as.data.frame()
          #rownames(annot) <-  var.mb[var]
          variant.ids.filtered.df.anno <- dplyr::bind_rows(variant.ids.filtered.df.anno, annot)  #readRDS('.//input/variant_ids_filtered_df_anno.rds')# 
        }
        
      
        # Original variant ids as rownames
        rownames(variant.ids.filtered.df.anno) <- variant.ids.filtered
      } else if(check_MBAPI() != 'MissionBio' && metadata()['genome_version',] == 'hg19'){
        # TODO
      } else if(metadata()['genome_version',] == 'hg38'){
        snpmart <- useEnsembl(biomart = "snp", dataset="hsapiens_snp")
        
        anno <- data.frame('refsnp_source' = c(), 'refsnp_id' = c(), 
                           'chr_name' = c(), 'chrom_start' = c(), 
                           'chrom_end' = c(), 'consequence_type_tv' = c(),
                           'clinical_significance' = c(), 'ensembl_gene_name' = c(), 'id' = c())
        
        for (var in variant.ids.filtered[1:3]){  
          message(paste0('annotate ', var))
          var.tmp <- str_match(var, "([0-9XY]+):(\\d+)") %>% 
            { paste(.[2], .[3], .[3], sep = ":") }
          tryCatch({
            anno.tmp <- getBM(
              attributes = c("refsnp_source", 'refsnp_id', 'chr_name', 'chrom_start', 'chrom_end',
                             "consequence_type_tv", "clinical_significance", 'ensembl_gene_name'),
              filters = 'chromosomal_region',
              values = var.tmp,
              mart = snpmart
            )
            anno.tmp$ID <- var
          }, error = function(e) {
            # Return NA for this SNP if an error occurs
            data.frame(refsnp_source = NA, refsnp_ID = NA, chr_name = NA, chrom_start = NA, 
                       chrom_end = NA, consequence_type_tv = NA, clinical_significance = NA, 
                       ensembl_gene_name = NA)
          })
          anno <- rbind(anno, anno.tmp)
        }
        
        lookup <- bitr(unique(anno$ensembl_gene_name), 'ENSEMBL', 'SYMBOL', org.Hs.eg.db)
        colnames(anno)[colnames(anno) == 'ensembl_gene_name'] <- 'ENSEMBL'
        variant.ids.filtered.df.anno <- merge(anno, lookup, all = T) %>%
          
          dplyr::select(any_of(c('ID', 'SYMBOL', 'Position', 'consequence_type_tv',
                                 'clinical_significance',
                                 'refsnp_id')))%>%
          dplyr::rename(any_of(c(
            Gene = 'SYMBOL',
            `chr_name` = 'Chr',
            `Consequence` = 'consequence_type_tv',
            `Clinical significance` = 'clinical_significance',
            dbsnp = 'refsnp_id'))) %>%
          tidyr::separate(ID, into = c('Chr', 'Start', 'Ref', 'Alt'), sep = ':', remove = F) %>% 
          as.data.frame()
        
      } else {
        message('Issue with MissionBio API. Using Biomart')
      }

      # Update variant ids transform to vector to get alphanumerical order of levels
      variant.ids.filtered.gene <- paste0(variant.ids.filtered.df.anno$Gene, ':', rownames(variant.ids.filtered.df.anno))
      variant.ids.filtered.gene <- factor(variant.ids.filtered.gene)
      colnames(vaf_matrix_filtered) <- variant.ids.filtered.gene
      colnames(genotype_matrix_filtered) <- variant.ids.filtered.gene
    })
    
    # Expose the visibility state
    output$plots_visible <- reactive({ plots_visible() })
    outputOptions(output, 'plots_visible', suspendWhenHidden=FALSE)
    
    
    # Render the filtered matrices and tables or plots as needed
    output$filtered_genotype <- renderTable({
      genotype_matrix_filtered
    })
    
    output$filtered_vaf <- renderTable({
      vaf_matrix_filtered  %>% datatable() %>%
        formatRound(columns = c(9, 12), digits = 2)
    })
    
    
    ## Variant panel plots -----------------------------------------------------
    ## Plot: No of variants ----------------------------------------------------
    output$var_plot1 <- renderPlot({
      req(plots_visible()) 
      
      plot(0,type='n',axes=FALSE,ann=FALSE)
      mtext(dim(vaf.matrix)[2], side = 3,line = -2, cex = 3, col = 'forestgreen')
      mtext('Number of variants total', side = 3, line = -4, cex = 1.5)
      
      # # Print ean mapped reads per cell
      mtext(dim(vaf_matrix_filtered)[2], side = 1,line = -4, cex = 3, col = 'dodgerblue')
      mtext('Number of variants filtered', side = 1, line = -2, cex = 1.5)
      
      # Draw box
      box(which = 'outer', lty = 'solid', col = 'grey')
    })
    
    ## Plot: No of cells ------------------------------------------------------
    output$var_plot2 <- renderPlot({
      req(plots_visible()) 
      
      # Print number of cells
      plot(0,type='n',axes=FALSE,ann=FALSE)
      mtext(dim(vaf.matrix)[1], side = 3,line = -2, cex = 3, col = 'forestgreen')
      mtext('Number of cells total', side = 3, line = -4, cex = 1.5)
      
      # # Print ean mapped reads per cell
      mtext(dim(vaf_matrix_filtered)[1], side = 1,line = -4, cex = 3, col = 'dodgerblue')
      mtext('Number of cells filtered', side = 1, line = -2, cex = 1.5)
      
      # Draw box
      box(which = 'outer', lty = 'solid', col = 'grey')
    })
    
    ## Heatmap: VAF I ---------------------------------------------------------
    output$var_plot3 <- renderPlot({
      colors.vaf <- circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292")) # TODO setup outside
     
      # Chromosome annotation    
      column_ha = HeatmapAnnotation(chr = factor(str_extract(colnames(vaf_matrix_filtered), 'chr(\\d|X|Y)+'),
                                                 levels = chromosomes),
                                    col = list(chr = chr_palette))
      
      vaf.matrix.filtered.hm <- vaf_matrix_filtered
      collect <- data.frame(row.names = '')
      
      # GT matrix annotation
      df <- do.call(rbind, lapply(genotype_matrix_filtered, function(x) {
        length(x) <- 4
        return(x)
      }))
      
      gt.anno <- data.frame(WT = integer(),
                            Het = integer(),
                            Hom = integer(),
                            Missing = integer())
      for (col in 1:ncol(genotype_matrix_filtered)){
        wt <- sum(genotype_matrix_filtered[,col] == 0)
        het <- sum(genotype_matrix_filtered[,col] == 1)
        hom <- sum(genotype_matrix_filtered[,col] == 2)
        mis <- sum(genotype_matrix_filtered[,col] == 3)
        gt.anno[col,] <- c(hom, het, wt,mis)
      }
      gt.anno$Total <- rowSums(gt.anno)

      proportions <- gt.anno %>%
        dplyr::mutate(across(c(Hom, Het, WT, Missing), ~ . / Total * 100)) %>%
        dplyr::select(-Total)
      
      #saveRDS(proportions, './input/prop_tmp.rds')
      
      rownames(proportions) <- rownames(variant.ids.filtered.df.anno)
      
      anno.bar <- anno_barplot(proportions, bar_width = 1, height = unit(3, "cm"),
                                gp = gpar(fill =  c("#D44292",
                                                    "#F6A97A",
                                                    "#414487FF",
                                                    "grey")))
      
      # Chromosome annotation and gt anno  
      column_ha <- HeatmapAnnotation(
        chr = factor(str_extract(colnames(vaf_matrix_filtered), 'chr(\\d|X|Y)+'),
                     levels = chromosomes),
        col = list(chr = chr_palette),
        'GT (%)' = anno.bar)
      
      colnames(vaf.matrix.filtered.hm) <- paste0(variant.ids.filtered.df.anno$Gene, ':', rownames(variant.ids.filtered.df.anno))
      vaf.matrix.filtered.hm <- vaf.matrix.filtered.hm
      draw(Heatmap(matrix = vaf.matrix.filtered.hm, 
                   name = 'VAF', 
                   col = colors.vaf,
                   show_column_dend = TRUE,
                   show_row_dend = FALSE, 
                   column_title = 'Filtered Variants',
                   row_title = 'Cells',  
                   top_annotation = column_ha,
           column_names_max_height=unit(100, "cm")))
    })
    
    output$legend <- renderPlot({
      # prep legend
      labels <- c("Hom", "Het", "WT", "Missing")
      colors <- c(
        `Hom` = "#D44292",
        `Het` = "#F6A97A",
        `WT` = "#414487FF",
        `Missing` = "#868686FF"
      )
      
      # Create the legend
      lgd <- Legend(labels = labels, 
                    title = "Genotype", 
                    legend_gp = gpar(fill = colors))
      draw(lgd)
    })
    
    ## Piechart: GT? ----------------------------------------------------------
    genotype_matrix_filtered <- genotype_matrix_filtered
    output$var_plot4 <- renderPlot({
      var_plot4 <- genotype_matrix_filtered %>%
        table() %>%
        as.data.frame() %>%
        rename(Genotype = '.') %>%  # Renaming the first column to 'Genotype'
        dplyr::mutate(
          Genotype = factor(dplyr::case_when(
            Genotype == 0 ~ 'WT',
            Genotype == '1' ~ 'Hom',
            Genotype == '2' ~ 'Het',
            TRUE ~ 'Missing'),
            levels = c('Hom', 'Het', 'WT', 'Missing'))) %>%
        mutate(prop = round((Freq / sum(Freq) * 100))) %>%
        dplyr::arrange(desc(Genotype)) %>%  # Sort Genotype in descending order for correct pie order
        mutate(cumulative = cumsum(prop), 
               lab.ypos = cumulative - 0.5 * prop) %>% 
        ggplot(aes(x = 2, y = prop, fill = Genotype)) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar(theta = "y", start = 0) +
        geom_text(aes(y = lab.ypos, label = paste0(prop, '%')), 
                  color = "white", size = 5) +
        xlim(0.5, 2.5)   +
        scale_fill_manual(values = mycols) +
        theme_default() +
        xlim(0.5, 2.5) +
        labs(title = 'Genotype Distribution (%)') +
        theme(axis.text = element_blank(),
              axis.title = element_blank())
      #ggsave(filename = './input/varplot_pie_dev.svg', device = grDevices::svg)
      var_plot4
    })
    
    ## Violin: GQ -------------------------------------------------------------
    output$var_plot5 <- renderPlot({
      colnames(genotype_matrix_filtered) <- rownames(variant.ids.filtered.df.anno)
      colnames(genoqual_matrix_filtered) <- rownames(variant.ids.filtered.df.anno)
      tmp.1 <- genotype_matrix_filtered %>% as.data.frame() %>% t() %>% melt(varnames = c('Variant', 'Cell'), value.name = 'Genotype')
      tmp.2 <-genoqual_matrix_filtered %>% as.data.frame() %>% t() %>% melt(varnames = c('Variant', 'Cell'), value.name = 'Genotype Quality')
      gt.df <- merge(tmp.1, tmp.2) %>%   
        dplyr::mutate(Genotype = ifelse(Genotype == 0, 'WT', 
                                        ifelse(Genotype == '1', 'Hom', 
                                               ifelse(Genotype == '2', 'Het', 'Missing'))))
      gt.df$Genotype <- factor(gt.df$Genotype, levels = c('Hom', 'Het', 'WT', 'Missing'))      
      
      var_plot5 <- ggplot(gt.df) +
        geom_violin(aes(x = Genotype, y =`Genotype Quality`, fill = Genotype), color = NA) +
        scale_fill_manual(values = mycols) +
        theme_default() +
        labs(title = 'Genotype Quality per Genotype (GATK)', x = NULL) +
        theme(panel.grid = element_blank()) +
        geom_hline(yintercept = 30, linetype = 'dashed')
      #ggsave(filename = './input/varplot_vio_dev.svg', device = grDevices::svg)
      var_plot5
      
    })
    
    # Variant tables ----------------------------------------------------------
    ## DT: overview filtered variants -----------------------------------------
    output$data_table_var <- renderDataTable({
      variant.ids.filtered.df.anno %>% 
        replace(is.na(.), '-') %>% 
        rownames_to_column(var = 'Variant') %>% 
        tidyr::separate(Variant, c('Chromosome', 'Position', 'Alt', 'Ref'), sep = ':|/', remove = F) %>%
        dplyr::mutate(Protein = ifelse(str_detect(Protein, '\\?'), Gene, Protein)) %>%
        dplyr::relocate(Gene, .before = Chromosome) %>%
        arrange(Gene) %>%
        datatable(.,  extensions = 'Buttons',
                  options = list(pageLength = 25, width = '95%',
                                 dom = 'Bfrtip',
                                 buttons = list( 
                                   list(extend = 'csv',   filename =  paste0("scafari_variants_", metadata()['sample_name',])),
                                   list(extend = 'excel', filename =  paste0("scafari_variants_", metadata()['sample_name',])),
                                   list(extend = 'pdf', filename =  paste0("scafari_variants_", metadata()['sample_name',])),
                                   list(extend = 'copy', filename =  paste0("scafari_variants_", metadata()['sample_name',])))),
                  rownames = F)
    })
    
    
    ## DT: select variants of interest ----------------------------------------
    output$data_table_var2 <- renderDataTable({

      sort(variant.ids.filtered.gene) %>%
        as.data.frame() %>%
        datatable(., 
                  options = list(pageLength = 25, width = '95%',
                                 scrollY = '600px',  
                                 paging = FALSE,      
                                 info = FALSE,
                                 dom = 'tf' 
                  ), rownames = F)
    })
    
    # Text: selected variants of interest --------------------------------------
    #? still included?
    output$selected_rows <- renderPrint({
      cat((sort(variant.ids.filtered.gene) %>% 
             as.data.frame() %>% 
             mutate(id = paste0(Gene, ':', rownames(.))) %>%
             dplyr::select(id) %>%
             arrange(id))[input$data_table_var2_rows_selected,])
    })
    
    # Plots explore panel ------------------------------------------------------
    ## Heatmap: VAF II ---------------------------------------------------------
    output$hm_1 <- renderPlot({
      # TODO before
      colors.vaf <- circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292"))  # TODO setup outside
      
      # Chromosome annotation    
      column_ha = HeatmapAnnotation(chr = factor(str_extract(colnames(vaf_matrix_filtered), 'chr(\\d|X|Y)+'),
                                                 levels = chromosomes),
                                    col = list(chr = chr_palette)
      )
      
      vaf.matrix.filtered.hm <- vaf_matrix_filtered  # TODO reienfolge prüfen
      collect <- data.frame(row.names = '')
      
      # GT matrix annotation
      df <- do.call(rbind, lapply(genotype_matrix_filtered, function(x) {
        # Ensure the vector is of length 4, filling with NA if necessary
        length(x) <- 4
        return(x)
      }))
      
      gt.anno <- data.frame(WT = integer(),
                            Het = integer(),
                            Hom = integer(),
                            Missing = integer())
      for (col in 1:ncol(genotype_matrix_filtered)){
        wt <- sum(genotype_matrix_filtered[,col] == 0)
        het <- sum(genotype_matrix_filtered[,col] == 1)
        hom <- sum(genotype_matrix_filtered[,col] == 2)
        mis <- sum(genotype_matrix_filtered[,col] == 3)
        gt.anno[col,] <- c(wt, het, hom, mis)
      }
      gt.anno$Total <- rowSums(gt.anno)
      proportions <- gt.anno %>%
        dplyr::mutate(across(c(WT, Het, Hom, Missing), ~ . / Total * 100)) %>%
        dplyr::select(-Total)# TODO check
      

      
      rownames(proportions) <- rownames(variant.ids.filtered.df.anno)
      
      anno.bar <- anno_barplot(proportions, bar_width = 1, height = unit(3, "cm"),
                                gp = gpar(fill =  c("#414487FF",
                                                    "#F6A97A",
                                                    "#D44292",
                                                    "grey")))
      
      # Chromosome annotation and gt anno  
      column_ha = HeatmapAnnotation(
        chr = factor(str_extract(colnames(vaf_matrix_filtered), 'chr(\\d|X|Y)+'),
                     levels = chromosomes),
        col = list(chr = chr_palette),
        'GT (%)' = anno.bar)

      colnames(vaf.matrix.filtered.hm) <- paste0(variant.ids.filtered.df.anno$Gene, ':', rownames(variant.ids.filtered.df.anno))
      draw(Heatmap(matrix = vaf.matrix.filtered.hm, 
                   name = 'VAF', 
                   col = colors.vaf,
                   show_column_dend = TRUE,
                   show_row_dend = FALSE, 
                   column_title = 'Filtered Variants',
                   row_title = 'Cells',  
                   top_annotation = column_ha,
                   column_names_max_height=unit(100, "cm")))
    })
    
    
    
    # Reactive variable to store user selections
    selected_variants <- reactiveVal(NULL)
    
    # Render the heatmap only when submit_var is clicked and selected_variants is not NULL
    # Update selected variants only when the button is clicked
    observeEvent(input$submit_var, {
      selected_variants(input$data_table_var2_rows_selected)  # Update the selected variants
    })
    
    # Render the heatmap only when submit_var is clicked and selected_variants is not NULL
    observeEvent(input$submit_var, {
      # Actualize heatmap
      output$hm_1 <- renderPlot({
        req(selected_variants())
        
        selected_variants_id <- order(variant.ids.filtered.gene)[selected_variants()]
        
        # Make colorpalette
        chromosomes <- c(paste0("chr", 1:21), "chrX", "chrY")
        colors.vaf <- circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292"))  # TODO outside
        
        vaf.matrix.filtered.hm <- vaf_matrix_filtered[, selected_variants_id]  # Directly use selected_variants
        
        column_ha = HeatmapAnnotation(chr = factor(str_extract(colnames(vaf.matrix.filtered.hm), 'chr(\\d|X|Y)+'),
                                                   levels = chromosomes),
                                      col = list(chr = chr_palette))
        
        collect <- data.frame(row.names = '')
        
        # GT matrix annotation (customize as needed)
        df <- do.call(rbind, lapply(genotype_matrix_filtered, function(x) {
          length(x) <- 4
          return(x)
        }))
        
        # Transform numerical genotype to WT, Het, Ho,, Missing dataframe
        gt.anno <- data.frame(WT = integer(), Het = integer(), Hom = integer(), Missing = integer())
        for (col in 1:ncol(genotype_matrix_filtered)){
          wt <- sum(genotype_matrix_filtered[, col] == 0)
          het <- sum(genotype_matrix_filtered[, col] == 1)
          hom <- sum(genotype_matrix_filtered[, col] == 2)
          mis <- sum(genotype_matrix_filtered[, col] == 3)
          gt.anno[col, ] <- c(wt, het, hom, mis)
        }
        
        # Force Het, Hom, Missing in the right order
        gt.anno$Total <- rowSums(gt.anno)
        proportions <- gt.anno %>%
          dplyr::mutate(across(c(WT, Het, Hom, Missing), ~ . / Total * 100)) %>%
          dplyr::select(-Total)
        rownames(proportions) <- rownames(variant.ids.filtered.df.anno)
        
        # Create and render the Heatmap
        Heatmap(matrix = vaf.matrix.filtered.hm, 
                name = 'VAF', 
                col = colors.vaf,
                show_column_dend = TRUE,
                show_row_dend = FALSE, 
                column_title = 'Filtered Variants',
                row_title = 'Cells',  
                top_annotation = column_ha,
                column_names_max_height=unit(100, "cm"))
      })
    })
    
    ## "Continue with selection" analysis --------------------------------------
    # Reactive variable for current variants
    current_variants <- reactiveVal(NULL)
    
    # Reactive variable for kneeplot data
    kneeplot_data <- reactiveVal(NULL)
    print.var <- reactiveVal(NULL)
    
    # Observer für den continue_var Button
    observeEvent(input$continue_var, {
      ## Elbow plot preparation ------------------------------------------------
      req(selected_variants()) 
      
      # Update current_variants with the selected variants
      current_variants(selected_variants())  
      
      current_variant_ids <-  sort(variant.ids.filtered.gene)[current_variants()]
      continue(TRUE)
      output$continue <- reactive({ continue() })
      outputOptions(output, 'continue', suspendWhenHidden=FALSE)
      
      # Compute the kneeplot data only upon clicking continue_var
      df <- vaf_matrix_filtered[, current_variant_ids] 
      df <- na.omit(df)
      df <- scale(df)
      
      # Determining Optimal Clusters - Elbow
      set.seed(123)
      plot <- fviz_nbclust(df, kmeans, method = "wss")
      
      # Store the plot in the reactive variable
      kneeplot_data(plot)
    })
    
    
    # Render the kneeplot using the reactive variable
    output$kneeplot <- renderPlot({
      req(kneeplot_data())  # Ensure there is data available
      # Render the kneeplot using the plot stored in kneeplot_data
      print(kneeplot_data())  # Render the plot
    })
    
    ## Clustering --------------------------------------------------------------
    # Define reactive variables for k2 and gg.clust
    k2 <- reactiveVal(NULL)
    gg.clust <- reactiveVal(NULL)
    ana_bar <- reactiveVal(NULL)
    cnv_plot4 <- reactiveVal(NULL)
    cnv_plot3 <- reactiveVal(NULL)
    cnv_plot6 <- reactiveVal(NULL)
    
    observe({
      if (is.numeric(input$n_clust) && input$n_clust >= 2) {
        enable("kmeans_btn")  # Enable button if valid
        runjs('document.getElementById("error_message").innerHTML = ""')  # Clear error message
      } else {
        disable("kmeans_btn")  # Disable button if not valid
        runjs('document.getElementById("error_message").innerHTML = "Please enter a numeric value greater than 2."')
      }
    })
    
    # Observe the k-means button event
    observeEvent(input$kmeans_btn, {
      req(current_variants())  # Ensure variants are selected
      req(is.numeric(input$n_clust) && input$n_clust > 2)  # Re-check the condition
      current_variant_ids <-  sort(variant.ids.filtered.gene)[current_variants()]

      # Print selected variants
      print.var(paste0("<ul>",
                       paste0(
                         "<li>",
                         current_variant_ids,
                         "</li>", 
                         collapse = ""
                       ),
                       "</ul>"))
      
      plots_visible_2(TRUE)
      
      
      ### Cluster plot ---------------------------------------------------------
      req(plots_visible_2())  # Ensure plots are visible and the data is available
      
      # Prepare the data
      df <- vaf_matrix_filtered[, current_variant_ids] #selected_variants()] 
      df <- na.omit(df)
      df <- scale(df)
      
      # Determining Optimal Clusters
      message('creating cluster plot')
      kmeans_result <- kmeans(df, centers = input$n_clust, nstart = 25) 
      
      # Store k2 in the reactive variable
      k2(kmeans_result)  # Store the kmeans result in the reactive variable
      
      # Generate the cluster plot and store it in shared_data
      gg.clust(fviz_cluster(kmeans_result, data = df, ellipse.type = "norm", geom = 'point') + 
                 theme_default())
      
      message('ggclust created')
      
      # Print the cluster plot for debugging
      print(gg.clust())  # Display the cluster plot
      
      gg.clust() #%>% saveRDS(., './input/cluster_plot.rds')

      
      ## Clustered Heatmap  ----------------------------------------------------
      req(gg.clust())  # Ensure there is a cluster plot available
      req(k2())  # Ensure k2 has a value
      
      # Make colorpalette
      chromosomes <- c(paste0("chr", 1:21), "chrX", "chrY")
      colors.vaf <- circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292"))  # TODO outsie
      
      # Here you don't need to call current_variants() or selected_variants() directly
      vaf.matrix.filtered.hm <- vaf_matrix_filtered[, current_variant_ids]  # Use current_variants directly
      column_ha <- HeatmapAnnotation(
        chr = factor(str_extract(colnames(vaf.matrix.filtered.hm), 'chr(\\d|X|Y)+'),
                     levels = chromosomes),
        col = list(chr = chr_palette)      )
      
      collect <- data.frame(row.names = '')
      
      # GT matrix annotation
      df <- do.call(rbind, lapply(genotype_matrix_filtered[, current_variant_ids], function(x) {
        length(x) <- 4
        return(x)
      }))
      
      # TODO is this used
      gt.anno <- data.frame(WT = integer(),
                            Het = integer(),
                            Hom = integer(),
                            Missing = integer())
      for (col in 1:ncol(genotype_matrix_filtered)){
        wt <- sum(genotype_matrix_filtered[,col] == 0)
        het <- sum(genotype_matrix_filtered[,col] == 1)
        hom <- sum(genotype_matrix_filtered[,col] == 2)
        mis <- sum(genotype_matrix_filtered[,col] == 3)
        gt.anno[col,] <- c(wt, het, hom, mis)
      }
      
      gt.anno$Total <- rowSums(gt.anno)
      proportions <- gt.anno %>%
        dplyr::mutate(across(c(WT, Het, Hom, Missing), ~ . / Total * 100)) %>%
        dplyr::select(-Total)# TODO check
      rownames(proportions) <- variant.ids.filtered.gene #rownames(variant.ids.filtered.df.anno)
      
      # Genotype annotation
      anno.bar = anno_barplot(proportions, bar_width = 1, height = unit(3, "cm"),
                              gp = gpar(fill =  c(WT = "#414487FF", 
                                                  Het = "#F6A97A", 
                                                  Hom = "#D44292",
                                                  Missing = "grey")))
      # Add cluster annoation
      colors <- gg_color_hue(input$n_clust)
      color_palette <- setNames(colors, as.character(seq(1,input$n_clust)))
      row_annot <- rowAnnotation(cluster = as.factor(k2()$cluster),
                                 col = list(cluster = color_palette))
      
      # Heatmap(matrix = vaf.matrix.filtered.hm, 
      #         name = 'VAF', 
      #         col = colors.vaf,
      #         show_column_dend = TRUE,
      #         show_row_dend = FALSE, 
      #         column_title = 'Filtered Variants',
      #         row_title = 'Cells',  
      #         top_annotation = column_ha,
      #         left_annotation = row_annot,
      #         row_split = as.factor(k2()$cluster))
      cnv_plot3(Heatmap(matrix = vaf.matrix.filtered.hm, 
                        name = 'VAF', 
                        col = colors.vaf,
                        show_column_dend = TRUE,
                        show_row_dend = FALSE, 
                        column_title = 'Filtered Variants',
                        row_title = 'Cells',  
                        top_annotation = column_ha,
                        left_annotation = row_annot,
                        row_split = as.factor(k2()$cluster),
                        column_names_max_height=unit(100, "cm")
      ))
      
      ## Violin: Explore variants ----------------------------------------------
      req(k2())
      req(plots_visible_2)
      
      vaf.matrix.filtered <- as.data.frame(vaf_matrix_filtered[,current_variant_ids] )
      colnames(vaf.matrix.filtered) <- current_variant_ids 
      
      # add cluster information
      vaf.matrix.filtered.tmp <- vaf.matrix.filtered
      
      colnames(vaf.matrix.filtered.tmp) <-current_variant_ids 
      vaf.matrix.filtered.tmp$cluster <- paste0('c', gg.clust()$data$cluster)
      vaf.matrix.filtered.tmp <- vaf.matrix.filtered.tmp %>%
        tidyr::pivot_longer(
          cols = c(-cluster),
          names_to = "variable",
          values_to = "value"
        ) %>% as.data.frame()
      vaf.matrix.filtered.tmp$variable <- factor(vaf.matrix.filtered.tmp$variable, levels = sort(current_variant_ids))
      
      cnv_plot4(vaf.matrix.filtered.tmp %>%
                  ggplot() +
                  geom_violin(aes(x = cluster, y = value, fill = cluster), alpha = 0.5, col = NA) +
                  geom_jitter(aes(x = cluster, y = value, col = cluster), size = 1) +
                  #geom_boxplot(aes(x = cluster, y = value), outliers = F) +
                  
                  #theme_default() +
                  labs(y = 'VAF', x = 'cluster')+
                  facet_grid(~variable) +
                  theme(
                    panel.grid = element_blank(),
                    panel.background = element_blank(),
                    title = element_text(size = 20),
                    text = element_text(size = 18)
                  ))
      
      ## Bar: Explore variants -------------------------------------------------
      req(k2())  # Ensure k2 has a value
      
      # TODO input several times
      genotype.matrix.filtered <- genotype_matrix_filtered[,current_variant_ids] 
      colnames(genotype.matrix.filtered) <- current_variant_ids   # TODO important, too
      
      # add cluster information
      genotype.matrix.filtered.tmp <- as.data.frame(genotype.matrix.filtered)
      
      genotype.matrix.filtered.tmp$cluster <- paste0('c', gg.clust()$data$cluster)
      genotype.matrix.filtered.tmp <<- melt(genotype.matrix.filtered.tmp)
      
      
      gt <- genotype.matrix.filtered.tmp %>%
        melt() %>%
        mutate(variable = factor(variable, levels = sort(current_variant_ids))) %>% 
        mutate(value = as.factor(value)) %>%
        mutate(Genotype = factor(dplyr::case_when(
          value == 0 ~ "WT",
          value == 1 ~ "Het",
          value == 2 ~ "Hom",
          TRUE ~ "Missing"  # Fallback for unexpected values
        ), levels = c('Hom', 'Het', 'WT', 'Missing')))
      
      
      recode_function <- function(x) {
        ifelse(x == 0, "WT",
               ifelse(x == 1, "Het",
                      ifelse(x == 2, "Hom",
                             ifelse(x == 3, "Missing", x))))
      }
      # Plot barplot
      ana_bar(gt %>%
                ggplot() +
                geom_bar(aes(x = cluster, fill = Genotype), col = NA) +
                #geom_jitter(aes(x = cluster, y = value, col = cluster), size = 1) +
                #geom_boxplot(aes(x = cluster, y = value), outliers = F) +
                scale_fill_manual(values = mycols.ngt) +
                #theme_default() +
                #labs(y = 'VAF', x = variant.of.interest)+
                facet_grid(~variable) +
                theme(
                  panel.grid = element_blank(),
                  panel.background = element_blank(),
                  title = element_text(size = 20),
                  text = element_text(size = 18)))
      
      
      ## Map: Explore variants colored by VAF ----------------------------------
      req(k2())
      req(plots_visible_2)
      
      rownames(vaf.matrix.filtered) <- paste0('cell', rownames(vaf.matrix.filtered)) # TODO cell ids
      gg.clust <- gg.clust()
      
      
      # TODO später nach cellID!
      rownames(gg.clust$data) <- paste0('cell', rownames(gg.clust$data))
      
      merged <- merge(gg.clust$data, vaf.matrix.filtered, by = 0)
      col_names <- colnames(merged)
      gene_cols <- col_names[grepl("[^:]:chr[[:digit:]]+", col_names, perl = T)]
      
      merged.var <- merged %>%
        tidyr::pivot_longer(
          cols = all_of(gene_cols),
          names_to = "variant", 
          values_to = "VAF" 
        )
      merged.var$variant <- factor(merged.var$variant, levels = sort(current_variant_ids))
      
      # Define the color function with gradient and -1 as grey
      color_func <- function(value) {
        ifelse(value == -1, "#BEBEBE",  # Grey for -1
               circlize::colorRamp2(c(0, 50, 100), c("#414487FF", "#F6A97A", "#D44292"))(value))}
      
      # Create a sequence of breakpoints covering the data range you want
      breakpoints <- seq(-1, 100, length.out = 101)  # Adjust based on your data, now includes -1
      
      # Generate the corresponding colors
      colors <- color_func(breakpoints)
      
      # legend
      data <- data.frame(
        value = c(-1, seq(0, 100, by = 25)),
        label = c("Missing", "0", "25", "50", "75", "100")
      )
     
      cnv_plot6(
        ggplot(merged.var, aes(x = x, y = y, color = VAF)) +
          
          geom_polygon(
            aes(x = x, y = y, group = cluster, fill = cluster), 
            alpha = 0.3,
            stat = "ellipse", 
            type = "norm", 
            level = 0.95, 
            segments = 51,
            na.rm = FALSE
          ) +
          geom_point() +
          scale_color_gradientn(colors = colors, values = scales::rescale(breakpoints)) +
          geom_point() +
          facet_grid(~factor(variant, levels = sort(current_variant_ids))) +
          theme(
            panel.grid = element_blank(),
            panel.background = element_blank(),
            title = element_text(size = 20),
            text = element_text(size = 18)) +
          labs(x = gg.clust$labels$x, y = gg.clust$labels$y)
      )
      
    })
    
    # End of Explore variants Panel plots
    # Plot calls Explore variants ----------------------------------------------
    output$k2_output <- renderPrint({
      req(k2())  # Ensure k2 has a value
      return(k2())  # Render the k-means result stored in k2
    })
    
    output$some_plot_output <- renderPlot({
      req(gg.clust())  # Ensure gg.clust has a value
      #saveRDS(gg.clust(),'./input/map_clust.rds')
      print(gg.clust())  # Render the plot stored in gg.clust
      
    })
    
    # Expose the visibility state
    output$plots_visible_2 <- reactive({ plots_visible_2() })
    outputOptions(output, 'plots_visible_2', suspendWhenHidden=FALSE)
    
    # Render the cluster plot using the reactive variable
    output$cluster_plot <- renderPlot({
      req(gg.clust())  # Ensure there is a cluster plot available
      req(plots_visible_2())  # Ensure plots are visible
      
      # Render the cluster plot using the plot stored in gg.clust
      print(gg.clust())  # Render the plot
    })
    
    output$cnv_plot3 <- renderPlot({
      print(cnv_plot3())
    })
    
    output$cnv_plot4 <- renderPlot({
      print(cnv_plot4())
    })
    
    
    output$cnv_plot6 <- renderPlot({
      print(cnv_plot6())
    })
    
    output$ana_bar <- renderPlot({
      print(ana_bar())
    })
    
    output$selected_rows_2 <- renderUI({
      HTML(print.var())
    })
  })
  
  # Create a waiter object
  w <- Waiter$new(id = "file", html = spin_3(), color = "rgba(255,255,255,0.8)")
  
  
  output$text1 <- renderText({paste("You have selected", input$var)})
  
  # Plots ----------------------------------------------------------------------
  # Sequencing plot 1
  output$seq_plot1 <- renderPlot({
    req(input$upload, metadata())
    validate(need(input$upload, 'Please, select a file to start'))  
    
    plot(0,type='n',axes=FALSE,ann=FALSE)
    mtext(metadata()['n_cells',], side = 3,line = -2, cex = 3, col = '#22A484FF')
    mtext('Number of cells', side = 3, line = -4, cex = 1.5)
    
    mtext(metadata()['sample_name',], side = 1,line = -4, cex = 3, col = '#22A484FF')  # TODO dyn
    mtext('Sample name', side = 1, line = -0.5, cex = 1.5)
    
    # Draw box
    box(which = 'outer', lty = 'solid', col = 'grey')
  })
  
  # Sequencing plot 2
  output$seq_plot2 <- renderPlot({
    req(input$upload)
    
    plot(0,type='n',axes=FALSE,ann=FALSE)
    mtext(paste0(round(as.numeric(metadata()['avg_panel_uniformity',])*100, digits = 2), ''), side = 3,line = -2, cex = 3, col = '#22A484FF')
    mtext('Panel uniformity (%)', side = 3, line = -4, cex = 1.5)
    
    # Print ean mapped reads per cell  # TODO check bei anderen, ob richtig!!
    mtext(round((as.numeric(metadata()['n_read_pairs_mapped_to_cells',])/as.numeric(metadata()['n_read_pairs',])*100), digits = 2), side = 1,line = -4, cex = 3, col = '#22A484FF')
    mtext('Read pairs assigned to cells (%)', side = 1, line = -2, cex = 1.5)
    
    # Draw box
    box(which = 'outer', lty = 'solid', col = 'grey')
    
  })
  
  # Sequencing log-log plot
  output$seq_plot3 <- renderPlot({
    req(input$upload, read.counts.df())
    
    counts.per.cell <- read.counts.df() %>% rowSums()
    
    seq_plot3 <- counts.per.cell %>% as.data.frame() %>%
      rownames_to_column('barcode') %>%
      ggplot(data = .) +
      geom_point(aes(x = reorder(barcode, .),y  = .),stat = 'identity') +
      labs(x = 'barcodes', y = 'Number of reads', title = 'Reads per Barcode') +
      theme_default() +
      theme(axis.line = element_line(),
            axis.ticks.y = element_line(),
            panel.grid = element_blank(),
            axis.text.x = element_blank()) +
      scale_y_log10() +
      scale_x_discrete(limits=rev)

    
    seq_plot3
  })
  
  # Panel plot 1
  output$panel_plot1 <- renderPlot({
    req(input$upload)
    
    plot(0,type='n',axes=FALSE,ann=FALSE)
    mtext(metadata()['panel_name',], side = 3,line = -2, cex = 3, col = '#22A884FF')
    mtext('Panel used', side = 3, line = -4, cex = 1.5)
    
    
    # Print ean mapped reads per cell
    mtext(metadata()['n_amplicons',], side = 1,line = -4, cex = 3, col = '#22A884FF')
    mtext('Number of amplicons', side = 1, line = -2, cex = 1.5)
    box(which = 'outer', lty = 'solid', col = 'grey')
    
  })
  
  # Panel plot 2
  output$panel_plot2 <- renderPlot({
    req(input$upload)
    req(amplicons())
    genes <- sapply(str_split(amplicons(), "_"), function(x) x[3])
    plot(0,type='n',axes=FALSE,ann=FALSE)
    mtext(length(unique(genes)), side =3,line = -2, cex = 3, col = '#F66D7A')
    mtext('Number of Genes covered', side = 3, line = -4, cex = 1.5)
    
    # Print ean mapped reads per cell
    mtext(floor(as.numeric(metadata()['n_read_pairs_mapped_to_cells',])/
                  as.numeric(metadata()['n_cells',])/as.numeric(metadata()['n_amplicons',])), side = 1,line = -4, cex = 3, col = '#22A484FF')  # TODO dyn
    mtext('Average read pairs per \namplicon and cell', side = 1, line = -0.5, cex = 1.5)
    box(which = 'outer', lty = 'solid', col = 'grey')
    #https://support.missionbio.com/hc/en-us/articles/360044185393-What-causes-a-low-of-mapped-reads-to-cells
  })
  
  # Panel karyoplot
  output$panel_plot3 <- renderPlot({
    # Prepare amplicon data
    gene.anno.gr <- gene.anno.df() %>%
      dplyr::mutate(Gene = str_split_i(id, '_', 3)) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)
    
    kp <- plotKaryotype(plot.type=1, cex = 1.5)
    kpPlotRegions(kp, data = gene.anno.gr, col="#F66D7A", data.panel = 1)
    
    # To avoid overplotting just plot each gene name once
    amp.genes <- gene.anno.df()  %>%
      dplyr::mutate(Gene = str_split_i(id, '_', 3)) %>%
      group_by(Gene) %>%
      slice(1) %>%
      dplyr::ungroup() %>%
      makeGRangesFromDataFrame(keep.extra.columns = T)
    
    # Add gene labels
    kpText(kp, chr=seqnames(amp.genes), x=start(amp.genes), y=0.75, labels=amp.genes$Gene, cex = 0.95)
  })
  
  output$panel_plot4 <- renderPlot({
    # Normalized mean read counts per amplicon
    rc <- as.data.frame(read.counts.df.norm() %>% colMeans()) %>% rownames_to_column('Amplicon')
    colnames(rc)[2] <- 'Normalized mean read\ncounts per amplicon'
    font_size <- max(3, min(18, 3 + (nrow(rc) - 1) * (2 / (3000 - 1))))
    panel_plot4 <- rc %>% ggplot(.) +
      geom_bar(aes(x = reorder(Amplicon, -`Normalized mean read\ncounts per amplicon`), y = `Normalized mean read\ncounts per amplicon`), stat = 'identity') +
      labs(x = '') +
      theme_default() +
      theme(axis.text.x = element_text(angle = -270, size = 10),
            axis.line = element_line())
    panel_plot4
  })
  
  output$panel_plot5 <- renderPlot({
    # Color panel uniformity
    # Panel uniformity is the percentage of targets that meet at least 20% of the
    # average depth of coverage.
    # Order by mean of read counts
    req(read.counts.df.norm(), gene.anno.df(), gene.anno.df()$id)
    read.counts.df.norm <- read.counts.df.norm()
    colnames(read.counts.df.norm) <- gene.anno.df()$id
    read.counts.df.norm.tmp <- read.counts.df.norm %>% melt()
    tmp.m <- read.counts.df.norm() %>% colMeans()
    tmp.m.o <- order(tmp.m)
    
    read.counts.df.norm.tmp$variable <- factor(read.counts.df.norm.tmp$variable, levels = names(tmp.m[tmp.m.o]))
    
    # Get mean average depth of coverage
    mean.tmp <- mean(as.matrix(read.counts.df.norm))
    message(paste0('Mean coverage: ', mean.tmp))
    
    # Extract amplicons with a lower mean coverage than the mean
    low.uniformity.amps <- names(tmp.m[tmp.m < 0.2*mean.tmp])
    
    # Introduce new column with uniformity info for color
    read.counts.df.norm.tmp <- read.counts.df.norm.tmp %>% dplyr::mutate(coverage = ifelse(variable %in% low.uniformity.amps, F, T))
    
    ggplot(read.counts.df.norm.tmp, aes(x = value, y = variable, color = coverage)) +
      geom_point() +
      scale_color_manual(values = c('TRUE' = '#414487FF', 'FALSE' = '#F66D7A')) +
      labs(x = 'Normalized read counts', y = 'Amplicons') +
      scale_x_continuous(limits = c(0, max(read.counts.df.norm.tmp$value)), expand = c(0,0)) +
      theme_default() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line = element_line(),
            legend.position = 'none')
  })
  
  
  # Tables ---------------------------------------------------------------------
  output$data_table_sample <- renderDataTable({
    (metadata() %>% rownames_to_column())[c(16, 14, 6, 5,4, 15),] %>% 
      mutate(rowname = gsub('n_', 'Number of ', rowname)) %>% 
      mutate(rowname = gsub('_', ' ', rowname)) %>%
      mutate(rowname = str_to_title(rowname)) %>% 
      datatable(., extension = 'Buttons', 
                options = list(pageLength = 6,
                               dom = 'Bt',
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_sequencing_sample_", metadata()['sample_name',])),
                                 list(extend = 'excel', filename =  paste0("scafari_sequencing_sample_", metadata()['sample_name',])),
                                 list(extend = 'pdf', filename =  paste0("scafari_sequencing_sample_", metadata()['sample_name',])),
                                 list(extend = 'copy', filename =  paste0("scafari_sequencing_sample_", metadata()['sample_name',])))),
                colnames = rep("", ncol(.)), rownames = F)
  })
  
  output$data_table_sequencing <- renderDataTable({
    rbind('Total read pairs' =  c(paste((as.numeric(metadata()['n_read_pairs',])))),
          'Read pairs trimmed' =  c(paste((as.numeric(metadata()['n_read_pairs_trimmed',])))),
          "Read pairs with valid barcodes" = c(paste(round(as.numeric(metadata()['n_read_pairs_valid_cell_barcodes',]))))) %>%
      
      datatable(.,  extensions = 'Buttons',
                options = list(pageLength = 5, width = '100%',
                               dom = 'Bt',
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_sequencing_overview_", metadata()['sample_name',])),
                                 list(extend = 'excel', filename =  paste0("scafari_sequencing_overview_", metadata()['sample_name',])),
                                 list(extend = 'pdf', filename =  paste0("scafari_sequencing_overview_", metadata()['sample_name',])),
                                 list(extend = 'copy', filename =  paste0("scafari_sequencing_overview_", metadata()['sample_name',])))),
                colnames = rep("", ncol(.)))
  })
  
  output$data_table_mapping <- renderDataTable({
    rbind('Reads mapped to genome (%)' =  c(round((as.numeric(metadata()['n_reads_mapped',])/(as.numeric(metadata()['n_read_pairs',])*2)*100),digits = 2)),
          'Reads mapped to target (%)' = c(round((as.numeric(metadata()['n_reads_mapped_insert',])/(as.numeric(metadata()['n_read_pairs',])*2)*100),digits = 2)))%>%
      datatable(., 
                extensions = 'Buttons',
                options = list(pageLength = 5, 
                               width = '100%',
                               dom = 'Bt',
                               buttons = list( 
                                 list(extend = 'csv',   filename =  paste0("scafari_sequencing_mapping", metadata()['sample_name',])),
                                 list(extend = 'excel', filename =  paste0("scafari_sequencing_mapping", metadata()['sample_name',])),
                                 list(extend = 'pdf', filename =  paste0("scafari_sequencing_mapping", metadata()['sample_name',])),
                                 list(extend = 'copy', filename =  paste0("scafari_sequencing_mapping", metadata()['sample_name',])))),
                class = 'display',
                colnames = rep("", ncol(.)))
  })
  
  ## Occurence of genes in panel -----------------------------------------------
  output$data_table_overview <- renderDataTable({
    if (metadata()['genome_version',] == 'hg19'){
      # Prepare exon database -----------------------------------------------------
      #Read Biomart Exon information and format them
      exons <- read.delim('./input//mart_export_grch37_p13.txt', header = T, sep = '\t')
      colnames(exons)[11] <- 'seqnames'
      colnames(exons)[5] <- 'start'
      colnames(exons)[6] <- 'end'
      colnames(exons)[7] <- 'Exon'
      exons$seqnames <- paste0('chr', exons$seqnames)
      exons.gr <- makeGRangesFromDataFrame(exons, keep.extra.columns = T)
      
      # Extract canonical transcripts
      known.canon <- read.delim('./input/UCSC_hg19_knownCanonical_chrom.bed', header = F, col.names = c('seqnames', 'start', 'end', 'Transcript.stable.ID.version', 'gene'))
      known.canon$Transcript.stable.ID.version <- gsub('\\..*', '', known.canon$Transcript.stable.ID.version)
      exons.gr$Transcript.stable.ID.version <- gsub('\\..*', '', exons.gr$Transcript.stable.ID.version)
      exons.gr.clean <- exons.gr[exons.gr$Transcript.stable.ID.version %in% known.canon$Transcript.stable.ID.version,]
      #TODO dynamic for different genome versions
      
      
      # Change chr info to merge them
      gene.anno.df.tmp <- gene.anno.df()
      #gene.anno.df.tmp$seqnames <- gsub('chr', '', gene.anno.df()$seqnames)
      gene.anno.df.tmp <- gene.anno.df.tmp %>% dplyr::mutate(Gene = str_split_i(id, '_', 3))
      gene.anno.gr <- makeGRangesFromDataFrame(gene.anno.df.tmp, keep.extra.columns = T)
      rm(gene.anno.df.tmp)
      exons.gr.clean <- exons.gr[exons.gr$Transcript.stable.ID.version %in%
                                   known.canon$Transcript.stable.ID.version,]
      ov <- findOverlaps(gene.anno.gr, exons.gr.clean, )
      
      # Define transcript columns
      mcols(gene.anno.gr)$`Exon` <- '-'
      mcols(gene.anno.gr)$`Canonical Transcript ID` <- '-'
      
      gene.anno.gr[queryHits(ov)]$`Exon` <- exons.gr.clean[subjectHits(ov)]$Exon  # i think its wrong
      gene.anno.gr[queryHits(ov)]$`Canonical Transcript ID` <- exons.gr.clean[subjectHits(ov)]$Transcript.stable.ID.version  # i think its wrong
      
      # Format for datatable
      df <- gene.anno.gr %>% as.data.frame()
      df %>% dplyr::select(seqnames, start, end, width, Gene, Exon, Canonical.Transcript.ID) %>%
        `colnames<-`(c('Chromosome', 'Start', 'End', 'Amplicon length (bp)', 'Gene', 'Exon', 'Canonical Transcript ID')) %>%
        datatable(., rownames = F,  extensions = 'Buttons',
                  options = list(pageLength = 10, width = '100%',
                                 dom = 'Bfrtip', 
                                 buttons = list( 
                                   list(extend = 'csv',   filename =  paste0("scafari_panel_", metadata()['sample_name',])),
                                   list(extend = 'excel', filename =  paste0("scafari_panel_", metadata()['sample_name',])),
                                   list(extend = 'pdf', filename =  paste0("scafari_panel_", metadata()['sample_name',])),
                                   list(extend = 'copy', filename =  paste0("scafari_panel_", metadata()['sample_name',])))))
      
    } else if (metadata()['genome_version',] == 'hg38'){
      message('hg38')
      
      withProgress(message = 'MANE annotation', detail = 'This may take a while', value = 0, {
        # MANE
        
        incProgress(1/3, 'Reading in MANE')
        mane.raw <- read.delim('./input/MANE.GRCh38.v1.3.ensembl_genomic.gff', skip = 2, header = F, sep = '\t')
        
        incProgress(1/3, 'Processing MANE')
        mane <- mane.raw %>% mutate(Exon = str_extract(V9, '(?<=exon_number=)\\d+'), 
                                    Gene = str_extract(V9, '(?<=gene_name=)[^;]+'),  
                                    `Transcript ID` = str_extract(V9, '(?<=transcript_id=)[^;]+'))  %>% 
          filter(V3 == 'exon')
        colnames(mane) <- c('seqnames','source', 'feature','start', 'end','score','strand', 'frame','Atrribute', 'Exon', 'Gene', 'Transcript ID')
        
        incProgress(1/5, 'Annotating')
        mane.gr <- makeGRangesFromDataFrame(mane, keep.extra.columns = T, na.rm = T)
        
        # Define transcript columns
        gene.anno.gr <- gene.anno.gr
        mcols(gene.anno.gr)[["Exon"]] <- '-'
        mcols(gene.anno.gr)[["Transcript ID"]] <- '-'
        mcols(gene.anno.gr)[["Gene"]] <- '-'
        ov <- findOverlaps(gene.anno.gr, mane.gr)
        
        mcols(gene.anno.gr)[queryHits(ov),]["Exon"] <- mcols(mane.gr)[subjectHits(ov),]["Exon"]
        mcols(gene.anno.gr)[queryHits(ov),]["Transcript ID"] <- mcols(mane.gr)[subjectHits(ov),]["Transcript ID"]
        mcols(gene.anno.gr)[queryHits(ov),]["Gene"] <- mcols(mane.gr)[subjectHits(ov),]["Gene"]
        
        gene.anno.gr %>% as.data.frame() %>% 
          select(id, seqnames, start, end, width, Gene, Exon, `Transcript.ID`) %>% 
          `colnames<-`(c('Amplicon ID', 'Chromosome', 'Start', 'End', 'Amplicon length (bp)', 'Gene', 'Exon', 'Canonical Transcript ID')) %>% 
          datatable(.,  extensions = 'Buttons',
                    options = list(pageLength = 10, width = '100%',
                                   dom = 'Bfrtip', 
                                   buttons = list( 
                                     list(extend = 'csv',   filename =  paste0("scafari_panel_", metadata()['sample_name',])),
                                     list(extend = 'excel', filename =  paste0("scafari_panel_", metadata()['sample_name',])),
                                     list(extend = 'pdf', filename =  paste0("scafari_panel_", metadata()['sample_name',])),
                                     list(extend = 'copy', filename =  paste0("scafari_panel_", metadata()['sample_name',])))), rownames = F)
      })
    } else {
      message('No proper genome version')
    }
  })
  
  session$onSessionEnded(function() {
    H5Fclose(h5_file())  # Close all open HDF5 handles
  })
  
}

# Run the application
shinyApp(ui = ui, server = server)
