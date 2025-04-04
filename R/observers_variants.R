createVariantFilteringObserver <- function(input, output, session, plots_visible) {
  
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
      variant.ids.filtered <<- variant.ids.filtered[!is.na(variant.ids.filtered)]
      rownames(vaf_matrix_filtered) <- NULL
    })
    
    ## Annotate variants -------------------------------------------------------
    # Test if missionbio api is available
    withProgress(message = 'Annotate Variants', value = 0, {
      
      if (check_MBAPI() == 'MissionBio' && metadata()['genome_version',] == 'hg19'){
        # Bring variants to MissionBio APIs variant format
        var.mb <- apply(variant.ids.filtered, 1, function(x){str_match(x, "(chr[0-9XY]+):(\\d+):([ATGC]+)/([ATGC]+)") %>%
            { paste(.[2], .[3], .[4], .[5], sep = "-") }})
        message(length(var.mb))
        variant.ids.filtered.df.anno <- data.frame(Gene = character(), Protein = character(), `Coding impact` = character(),
                                                   Function = character(), DANN = character(), ClinVar = character() , dbsnp = character())
        
        for (var in seq_along(var.mb)) {
          incProgress(1/length(var.mb), detail = paste0('Annotate variant ', var))
          
          # Reach MissioBio API
          url <- paste0("https://api.missionbio.io/annotations/v1/variants?ids=", var.mb[var])
          
          # Parse URL to get variant annotation
          res = httr::GET(url)
          data = jsonlite::fromJSON(rawToChar(res$content))  # TODO vary
          print(length(unlist(data$annotations)))
          annot <- data$annotations %>%
            unlist() %>%
            t() %>%
            as.data.frame() %>%
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
          print(var)
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
          print(anno.tmp)
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
      message()
      
      # Update variant ids transform to vector to get alphanumerical order of levels
      variant.ids.filtered.gene <- paste0(variant.ids.filtered.df.anno$Gene, ':', rownames(variant.ids.filtered.df.anno))
      variant.ids.filtered.gene <<- factor(variant.ids.filtered.gene)
      message('to factors done')
      colnames(vaf_matrix_filtered) <- variant.ids.filtered.gene
      message('colnames vafs done')
      colnames(genotype_matrix_filtered) <- variant.ids.filtered.gene
      message('colnames gt done')
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
        print(col)
        wt <- sum(genotype_matrix_filtered[,col] == 0)
        het <- sum(genotype_matrix_filtered[,col] == 1)
        hom <- sum(genotype_matrix_filtered[,col] == 2)
        mis <- sum(genotype_matrix_filtered[,col] == 3)
        gt.anno[col,] <- c(hom, het, wt,mis)
      }
      gt.anno$Total <- rowSums(gt.anno)
      gt.anno <- gt.anno
      
      proportions <- gt.anno %>%
        dplyr::mutate(across(c(Hom, Het, WT, Missing), ~ . / Total * 100)) %>%
        dplyr::select(-Total)
      
      saveRDS(proportions, './input/prop_tmp.rds')
      message(dim(variant.ids.filtered.df.anno))
      message(dim(proportions))
      
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
      message(paste0('chr_pal ', chr_palette))
      
      colnames(vaf.matrix.filtered.hm) <- paste0(variant.ids.filtered.df.anno$Gene, ':', rownames(variant.ids.filtered.df.anno))
      vaf.matrix.filtered.hm <<- vaf.matrix.filtered.hm
      draw(Heatmap(matrix = vaf.matrix.filtered.hm, 
                   name = 'VAF', 
                   col = colors.vaf,
                   show_column_dend = TRUE,
                   show_row_dend = FALSE, 
                   column_title = 'Filtered Variants',
                   row_title = 'Cells',  
                   top_annotation = column_ha))
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
      print(head(genotype_matrix_filtered))
      print(dim(genotype_matrix_filtered))
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
      ggsave(filename = './input/varplot_pie_dev.svg', device = grDevices::svg)
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
      ggsave(filename = './input/varplot_vio_dev.svg', device = grDevices::svg)
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
      message(paste0('!!!!!\nIDS: ', variant.ids.filtered.gene))
      
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
      message(input$data_table_var2_rows_selected)
      message(variant.ids.filtered.gene)
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
      
      vaf.matrix.filtered.hm <<- vaf_matrix_filtered  # TODO reienfolge prüfen
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
        print(col)
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
      
      
      message(dim(variant.ids.filtered.df.anno))
      message(dim(proportions))
      
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
      message(paste0('chr_pal ', chr_palette))
      
      colnames(vaf.matrix.filtered.hm) <- paste0(variant.ids.filtered.df.anno$Gene, ':', rownames(variant.ids.filtered.df.anno))
      draw(Heatmap(matrix = vaf.matrix.filtered.hm, 
                   name = 'VAF', 
                   col = colors.vaf,
                   show_column_dend = TRUE,
                   show_row_dend = FALSE, 
                   column_title = 'Filtered Variants',
                   row_title = 'Cells',  
                   top_annotation = column_ha))
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
                top_annotation = column_ha)
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
    
    # Observer für den k-means Button
    observeEvent(input$kmeans_btn, {
      req(current_variants())  # Ensure that there are selected variants
      
      current_variant_ids <-  sort(variant.ids.filtered.gene)[current_variants()]
      print(paste0('curr ids: ', current_variant_ids))
      
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
      message('ggclust printed')
      
      
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
        print(col)
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
      saveRDS(proportions, './input/debug_proportions.rds')
      
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
                        row_split = as.factor(k2()$cluster)
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
      genotype.matrix.filtered.tmp <<- as.data.frame(genotype.matrix.filtered)
      
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
}