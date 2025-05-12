# Server -----------------------------------------------------------------------
app_server <- function(input, output, session) {
  # Setup reactivity -----------------------------------------------------------
  plots_visible <- reactiveVal(FALSE)
  plots_visible_2 <- reactiveVal(FALSE)
  continue <- reactiveVal(FALSE)

  shinyjs::disable("filter_btn")
  shinyjs::disable("continue_var")
  shinyjs::disable("kmeans_btn")

  # Upload settings
  options(shiny.maxRequestSize = 500 * 1024^2)

  # End of setting up  ---------------------------------------------------------

  # Process input --------------------------------------------------------------
  ## Check input ----------------------------------------------------------------
  h5_file <- reactive({
    req(input$upload) # Ensure input$upload exists

    # Open the HDF5 file
    h5_file <- H5Fopen(input$upload$datapath)

    # Check existence of important slots in input
    check_result <- checkH5(h5_file)

    # If input is valid upload gets deactivated and buttons are activated
    if (check_result$valid) {
      shinyjs::enable("filter_btn") # Enable button if valid
      shinyjs::disable("upload") # Deactivate input if valid
      output$status <- renderText("File uploaded successfully!")
      return(h5_file)
    } else {
      showNotification(
        paste(
          "Error: The following paths are missing in the H5 file:",
          paste(check_result$missing, collapse = ", ")
        ),
        type = "error",
        duration = NULL
      )
      return(NULL) # Return NULL if invalid
    }
  })

  # Display conditional panels if file is uploaded successfully
  output$file_ready <- reactive({
    !is.null(h5_file())
  })
  outputOptions(output, "file_ready", suspendWhenHidden = FALSE)

  # Panel analyses--------------------------------------------------------------
  sce <- reactive({
    req(input$upload$datapath)
    sce_amp <- h5ToSce(input$upload$datapath)$sce_amp
    norm <- normalizeReadCounts(sce = sce_amp)
    return(norm)
  })

  se.var <- reactive({
    req(input$upload$datapath)
    return(h5ToSce(input$upload$datapath)$se_var)
  })

  # Variant analyses -----------------------------------------------------------
  ## Variant filtering ---------------------------------------------------------
  # Filter variants after filter button hit
  rv <- reactiveValues(sce_filtered = NULL)
  observeEvent(input$filter_btn, {
    plots_visible(TRUE)
    current_sce <- sce()
    se.var <- se.var()
    filteres <- filterVariants(
      depth.threshold = 10,
      genotype.quality.threshold = 30,
      vaf.ref = 5,
      vaf.het = 35,
      vaf.hom = 95,
      min.cell = 50,
      min.mut.cell = 1,
      se.var = se.var,
      sce = current_sce,
      shiny = TRUE
    )

    indices_to_keep <- match(filteres$cells.keep, colData(current_sce)[[1]],
      nomatch = 0
    )
    if (length(indices_to_keep) > 0) {
      se.f <- SummarizedExperiment(
        assays = list(
          VAF = t(filteres$vaf.matrix.filtered),
          Genotype = t(filteres$genotype.matrix.filtered),
          Genoqual = t(filteres$genoqual.matrix.filtered)
        ),
        rowData = filteres$variant.ids.filtered,
        colData = filteres$cells.keep
      )

      sce_filtered <- current_sce[, indices_to_keep, drop = FALSE]
      SingleCellExperiment::altExp(sce_filtered, "variants") <- se.f
      sce_filtered <- annotateVariants(sce = sce_filtered, shiny = TRUE)

      # Update the reactive value
      rv$sce_filtered <- sce_filtered
    } else {
      showModal(modalDialog(
        title = "Warning",
        "No cells passed the filtering criteria."
      ))
    }

    # Expose the visibility state
    output$plots_visible <- reactive({
      plots_visible()
    })
    outputOptions(output, "plots_visible", suspendWhenHidden = FALSE)

    ## Variant panel plots -----------------------------------------------------
    ## Plot: No of variants ----------------------------------------------------
    output$var_plot1 <- renderPlot({
      req(plots_visible())
      req(rv$sce_filtered)

      plot(0, type = "n", axes = FALSE, ann = FALSE)
      mtext(dim(se.var())[1], side = 3, line = -2, cex = 3, col = "forestgreen")
      mtext("Number of variants total", side = 3, line = -4, cex = 1.5)

      # # Print ean mapped reads per cell
      mtext(dim(altExp(rv$sce_filtered))[1],
        side = 1, line = -4, cex = 3,
        col = "dodgerblue"
      )
      mtext("Number of variants filtered", side = 1, line = -2, cex = 1.5)

      # Draw box
      box(which = "outer", lty = "solid", col = "grey")
    })

    ## Plot: No of cells ------------------------------------------------------
    output$var_plot2 <- renderPlot({
      req(plots_visible())

      # Print number of cells
      plot(0, type = "n", axes = FALSE, ann = FALSE)
      mtext(dim(sce())[2], side = 3, line = -2, cex = 3, col = "forestgreen")
      mtext("Number of cells total", side = 3, line = -4, cex = 1.5)

      # Print ean mapped reads per cell
      mtext(dim(altExp(rv$sce_filtered))[2],
        side = 1, line = -4, cex = 3,
        col = "dodgerblue"
      )
      mtext("Number of cells filtered", side = 1, line = -2, cex = 1.5)

      # Draw box
      box(which = "outer", lty = "solid", col = "grey")
    })

    ## Heatmap: VAF I ---------------------------------------------------------
    output$var_plot3 <- renderPlot({
      req(rv$sce_filtered)
      plotVariantHeatmap(rv$sce_filtered)
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
      lgd <- Legend(
        labels = labels,
        title = "Genotype",
        legend_gp = gpar(fill = colors)
      )
      draw(lgd)
    })


    ## Violin: GQ -------------------------------------------------------------
    output$var_plot5 <- renderPlot({
      plotGenotypequalityPerGenotype(sce_filtered) +
        theme(
          title = element_text(size = 20),
          text = element_text(size = 16)
        )
    })

    # Variant tables ----------------------------------------------------------
    ## DT: overview filtered variants -----------------------------------------
    output$data_table_var <- renderDataTable({
      sample.name <- sce_filtered@metadata[["sample_name"]]
      file.out <- paste0("scafari_variants_", sample.name)
      rowData(altExp(sce_filtered)) %>%
        as.data.frame() %>%
        dplyr::select(-any_of("id")) %>%
        replace(is.na(.), "-") %>%
        rownames_to_column(var = "Variant") %>%
        tidyr::separate(Variant, c("Chromosome", "Position", "Alt", "Ref"),
          sep = ":|/", remove = FALSE
        ) %>%
        dplyr::mutate(Protein = ifelse(str_detect(Protein, "\\?"),
          Gene, Protein
        )) %>%
        dplyr::relocate(Gene, .before = Chromosome) %>%
        arrange(Gene) %>%
        datatable(.,
          extensions = "Buttons",
          options = list(
            pageLength = 25, width = "95%",
            dom = "Bfrtip",
            buttons = list(
              list(extend = "csv", filename = file.out),
              list(extend = "excel", filename = file.out),
              list(extend = "pdf", filename = file.out),
              list(extend = "copy", filename = file.out)
            )
          ),
          rownames = FALSE
        )
    })


    ## DT: select variants of interest -----------------------------------------
    output$data_table_var2 <- renderDataTable({
      sort(paste0(
        rowData(altExp(sce_filtered))$Gene, ":",
        rowData(altExp(sce_filtered))$id
      )) %>%
        as.data.frame() %>%
        datatable(.,
          options = list(
            pageLength = 25, width = "95%",
            scrollY = "600px",
            paging = FALSE,
            info = FALSE,
            dom = "tf"
          ), rownames = FALSE
        )
    })

    # Text: selected variants of interest --------------------------------------
    output$selected_rows <- renderPrint({
      cat((sort(paste0(
        rowData(altExp(sce_filtered))$Gene, ":",
        rowData(altExp(sce_filtered))$id
      )) %>%
        as.data.frame() %>%
        mutate(id = paste0(Gene, ":", rownames(.))) %>%
        dplyr::select(id) %>%
        arrange(id))[input$data_table_var2_rows_selected, ])
    })

    # Plots explore panel ------------------------------------------------------
    ## Heatmap: VAF II ---------------------------------------------------------
    output$hm_1 <- renderPlot({
      req(rv$sce_filtered) # Ensures that the plot only renders when sce_filtered is available
      plotVariantHeatmap(rv$sce_filtered)
    })

    # Reactive variable to store user selections
    selected_variants <- reactiveVal(NULL)

    # Render the heatmap only when submit_var is clicked and selected_variants is not NULL
    # Update selected variants only when the button is clicked
    observeEvent(input$submit_var, {
      shinyjs::enable("continue_var")

      selected_variants(input$data_table_var2_rows_selected) # Update the selected variants
    })

    # Render the heatmap only when submit_var is clicked and selected_variants is not NULL
    observeEvent(input$submit_var, {
      # Actualize heatmap
      output$hm_1 <- renderPlot({
        req(selected_variants())
        req(rv$sce_filtered) # Ensures that the plot only renders when sce_filtered is available
        sce_filtered <- rv$sce_filtered
        output$hm_1 <- renderPlot({
          variant.ids.filtered.gene <- paste0(rowData(altExp(sce_filtered))$Gene, ":", rowData(altExp(sce_filtered))$id)
          selected_variants_id <- order(variant.ids.filtered.gene)[selected_variants()]

          # Load matrices
          vaf_matrix_filtered <- as.data.frame(t(assay(altExp(sce_filtered, "variants"), "VAF")))
          colnames(vaf_matrix_filtered) <- paste0(rowData(altExp(sce_filtered, "variants"))$Gene, ":", rowData(altExp(sce_filtered, "variants"))$id)
          genotype.matrix.filtered <- as.data.frame(t(assay(altExp(sce_filtered), "Genotype")))
          colnames(genotype.matrix.filtered) <- paste0(rowData(altExp(sce_filtered))$Gene, ":", rowData(altExp(sce_filtered))$id)
          row_data <- rowData(altExp(sce_filtered, "variants"))
          names <- paste0(row_data$Gene, ":", as.character(row_data$id))

          vaf.matrix.filtered.hm <- vaf_matrix_filtered[, selected_variants_id] # Directly use selected_variants

          column_ha <- HeatmapAnnotation(
            chr = factor(str_extract(colnames(vaf.matrix.filtered.hm), "chr(\\d|X|Y)+"),
              levels = chromosomes
            ),
            col = list(chr = chr_palette)
          )

          collect <- data.frame(row.names = "")

          # GT matrix annotation (customize as needed)
          df <- do.call(rbind, lapply(genotype.matrix.filtered, function(x) {
            length(x) <- 4
            return(x)
          }))

          # Transform numerical genotype to WT, Het, Ho,, Missing dataframe
          gt.anno <- data.frame(
            WT = integer(), Het = integer(),
            Hom = integer(), Missing = integer()
          )
          for (col in 1:ncol(genotype.matrix.filtered)) {
            wt <- sum(genotype.matrix.filtered[, col] == 0)
            het <- sum(genotype.matrix.filtered[, col] == 1)
            hom <- sum(genotype.matrix.filtered[, col] == 2)
            mis <- sum(genotype.matrix.filtered[, col] == 3)
            gt.anno[col, ] <- c(wt, het, hom, mis)
          }

          # Force Het, Hom, Missing in the right order
          gt.anno$Total <- rowSums(gt.anno)
          proportions <- gt.anno %>%
            dplyr::mutate(across(c(WT, Het, Hom, Missing), ~ . / Total * 100)) %>%
            dplyr::select(-Total)
          rownames(proportions) <- names #+rownames(variant.ids.filtered.df.anno)

          colors.vaf <- circlize::colorRamp2(
            c(0, 50, 100),
            c("#414487FF", "#F6A97A", "#D44292")
          )

          # Create and render the Heatmap
          Heatmap(
            matrix = vaf.matrix.filtered.hm,
            name = "VAF",
            col = colors.vaf,
            show_column_dend = TRUE,
            show_row_dend = FALSE,
            column_title = "Filtered Variants",
            row_title = "Cells",
            top_annotation = column_ha
          )
        })
      })
    })

    ## "Continue with selection" analysis --------------------------------------
    # Reactive variable for current variants
    current_variants <- reactiveVal(NULL)
    kneeplot_data <- reactiveVal(NULL)
    print.var <- reactiveVal(NULL)

    # Observer für den continue_var Button
    observeEvent(input$continue_var, {
      ## Elbow plot preparation ------------------------------------------------
      req(selected_variants())
      # TODO maybe define reactive
      variant.ids.filtered.gene <- paste0(
        rowData(altExp(sce_filtered))$Gene,
        ":", rowData(altExp(sce_filtered))$id
      )

      # Update current_variants with the selected variants
      current_variants(selected_variants())
      current_variant_ids <- sort(variant.ids.filtered.gene)[current_variants()]

      continue(TRUE)
      output$continue <- reactive({
        continue()
      })
      outputOptions(output, "continue", suspendWhenHidden = FALSE)

      # Compute the kneeplot data only upon clicking continue_var
      # Store the plot in the reactive variable
      kneeplot_data(plotElbow(sce_filtered, current_variant_ids))
    })


    # Render the kneeplot using the reactive variable
    output$kneeplot <- renderPlot({
      req(kneeplot_data()) # Ensure there is data available
      print(kneeplot_data()) # Render the plot
    })

    ## Clustering --------------------------------------------------------------
    # Define reactive variables for k2 and gg.clust
    k2 <- reactiveVal(NULL)
    gg.clust <- reactiveVal(NULL)
    ana_bar <- reactiveVal(NULL)
    vaf_hm <- reactiveVal(NULL)
    vaf_violin <- reactiveVal(NULL)
    vaf_map <- reactiveVal(NULL)

    observe({
      if (is.numeric(input$n_clust) && input$n_clust >= 2) {
        enable("kmeans_btn") # Enable button if valid
        runjs('document.getElementById("error_message").innerHTML = ""')
      } else {
        disable("kmeans_btn") # Disable button if not valid
        runjs('document.getElementById("error_message").innerHTML = "Please enter a numeric value greater than 2."')
      }
    })

    # Observe the k-means button event
    observeEvent(input$kmeans_btn, {
      req(current_variants()) # Ensure variants are selected
      req(is.numeric(input$n_clust) && input$n_clust >= 2) # Re-check the condition
      variant.ids.filtered.gene <- paste0(
        rowData(altExp(sce_filtered))$Gene,
        ":", rowData(altExp(sce_filtered))$id
      )
      variants.of.interest <- sort(variant.ids.filtered.gene)[current_variants()]

      # Print selected variants
      print.var(paste0(
        "<ul>",
        paste0(
          "<li>",
          variants.of.interest,
          "</li>",
          collapse = ""
        ),
        "</ul>"
      ))

      plots_visible_2(TRUE)

      ### Cluster plot ---------------------------------------------------------
      req(plots_visible_2()) # Ensure plots are visible and the data available
      cluster.res <- clusterVariantSelection(sce_filtered,
                                             variants.of.interest, 
                                             input$n_clust)
      k2(cluster.res[["k_means"]])
      gg.clust(cluster.res[["clusterplot"]])

      ## Clustered Heatmap  ----------------------------------------------------
      req(gg.clust()) # Ensure there is a cluster plot available
      req(k2())

      # Make colorpalette
      chromosomes <- c(paste0("chr", 1:21), "chrX", "chrY")
      colors.vaf <- circlize::colorRamp2(c(0, 50, 100), 
                                         c("#414487FF", "#F6A97A", "#D44292"))

      vaf.matrix.filtered <- as.data.frame(t(assay(altExp(sce_filtered, 
                                                          "variants"), "VAF")))
      colnames(vaf.matrix.filtered) <- 
        paste0(rowData(altExp(sce_filtered,"variants"))$Gene, 
               ":", rowData(altExp(sce_filtered, 
               "variants"))$id)

      genotype.matrix.filtered <- as.data.frame(t(assay(altExp(sce_filtered), 
                                                        "Genotype")))
      colnames(genotype.matrix.filtered) <- 
        paste0(rowData(altExp(sce_filtered))$Gene, ":", 
               rowData(altExp(sce_filtered))$id)

      vaf.matrix.filtered.hm <- vaf.matrix.filtered[, variants.of.interest]
      column_ha <- HeatmapAnnotation(
        chr = factor(str_extract(colnames(vaf.matrix.filtered.hm), 
                                 "chr(\\d|X|Y)+"),
          levels = chromosomes),
        col = list(chr = chr_palette))
      collect <- data.frame(row.names = "")

      # GT matrix annotation
      df <- do.call(rbind, lapply(genotype.matrix.filtered[, variants.of.interest], function(x) {
        length(x) <- 4
        return(x)
      }))

      gt.anno <- data.frame(
        WT = integer(),
        Het = integer(),
        Hom = integer(),
        Missing = integer()
      )
      
      for (col in 1:ncol(genotype.matrix.filtered)) {
        wt <- sum(genotype.matrix.filtered[, col] == 0)
        het <- sum(genotype.matrix.filtered[, col] == 1)
        hom <- sum(genotype.matrix.filtered[, col] == 2)
        mis <- sum(genotype.matrix.filtered[, col] == 3)
        gt.anno[col, ] <- c(wt, het, hom, mis)
      }

      gt.anno$Total <- rowSums(gt.anno)
      proportions <- gt.anno %>%
        dplyr::mutate(across(c(WT, Het, Hom, Missing), ~ . / Total * 100)) %>%
        dplyr::select(-Total)
      rownames(proportions) <- variant.ids.filtered.gene

      # Genotype annotation
      anno.bar <- anno_barplot(proportions,
        bar_width = 1, height = unit(3, "cm"),
        gp = gpar(fill = c(
          WT = "#414487FF",
          Het = "#F6A97A",
          Hom = "#D44292",
          Missing = "grey"
        ))
      )
      # Add cluster annoation
      colors <- gg_color_hue(input$n_clust)
      color_palette <- setNames(colors, as.character(seq(1, input$n_clust)))
      row_annot <- rowAnnotation(
        cluster = as.factor(k2()$cluster),
        col = list(cluster = color_palette)
      )
      vaf_hm(Heatmap(
        matrix = vaf.matrix.filtered.hm,
        name = "VAF",
        col = colors.vaf,
        show_column_dend = TRUE,
        show_row_dend = FALSE,
        column_title = "Filtered Variants",
        row_title = "Cells",
        top_annotation = column_ha,
        left_annotation = row_annot,
        row_split = as.factor(k2()$cluster)
      ))

      ## Violin: Explore variants ----------------------------------------------
      req(k2())
      req(plots_visible_2)
      violin <- plotClusterVAF(sce_filtered,
        variants.of.interest = variants.of.interest,
        gg.clust = gg.clust()
      ) +
        theme(
          title = element_text(size = 20),
          text = element_text(size = 16)
        )
      vaf_violin(violin)

      ## Bar: Explore variants -------------------------------------------------
      req(k2())
      ana_bar(plotClusterGenotype(sce_filtered,
        variants.of.interest = variants.of.interest,
        gg.clust = gg.clust()
      ) +
        theme(
          title = element_text(size = 20),
          text = element_text(size = 16)
        ))

      ## Map: Explore variants colored by VAF ----------------------------------
      req(k2())
      sce_filtered <- sce_filtered
      variants.of.interest <- variants.of.interest
      gg.clust <- gg.clust()
      vaf_map(plotClusterVAFMap(sce_filtered,
        variants.of.interest = variants.of.interest,
        gg.clust = gg.clust()
      ) +
        theme(
          title = element_text(size = 20),
          text = element_text(size = 16)
        ))
    })

    # End of Explore variants Panel plots
    # Plot calls Explore variants ----------------------------------------------
    output$k2_output <- renderPrint({
      req(k2())
      return(k2()) # Render the k-means result stored in k2
    })

    # Expose the visibility state
    output$plots_visible_2 <- reactive({
      plots_visible_2()
    })
    outputOptions(output, "plots_visible_2", suspendWhenHidden = FALSE)

    # Render the cluster plot using the reactive variable
    output$cluster_plot <- renderPlot({
      print(gg.clust())
    })

    output$vaf_hm <- renderPlot({
      print(vaf_hm())
    })

    output$vaf_violin <- renderPlot({
      print(vaf_violin())
    })

    output$vaf_map <- renderPlot({
      print(vaf_map())
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

  output$text1 <- renderText({
    paste("You have selected", input$var)
  })


  # Plots ----------------------------------------------------------------------
  # Sequencing plot 1
  output$seq_plot1 <- renderPlot({
    req(input$upload, sce())
    sce_obj <- sce()
    validate(need(input$upload, "Please, select a file to start"))

    metadata <- sce_obj@metadata

    plot(0, type = "n", axes = FALSE, ann = FALSE)
    sample_name <- metadata[["sample_name"]]
    n_cells <- metadata[["n_cells"]]

    mtext(n_cells, side = 3, line = -2, cex = 3, col = "#22A884FF")
    mtext("Number of cells", side = 3, line = -4, cex = 1.5)

    mtext(sample_name, side = 1, line = -4, cex = 3, col = "#22A884FF")
    mtext("Sample name", side = 1, line = -0.5, cex = 1.5)

    # Draw box
    box(which = "outer", lty = "solid", col = "grey")
  })

  # Sequencing plot 2
  output$seq_plot2 <- renderPlot({
    req(input$upload)
    sce_obj <- sce()
    validate(need(input$upload, "Please, select a file to start"))
    metadata <- sce_obj@metadata

    plot(0, type = "n", axes = FALSE, ann = FALSE)
    mtext(
      paste0(round(as.numeric(metadata[["avg_panel_uniformity"]]) * 100,
        digits = 2
      ), ""),
      side = 3, line = -2, cex = 3,
      col = "#22A484FF"
    )
    mtext("Panel uniformity (%)", side = 3, line = -4, cex = 1.5)

    # Print ean mapped reads per cell  # TODO check bei anderen, ob richtig!!
    mtext(
      round((as.numeric(metadata[["n_read_pairs_mapped_to_cells"]]) /
        as.numeric(metadata[["n_read_pairs"]]) * 100), digits = 2),
      side = 1, line = -4, cex = 3, col = "#22A484FF"
    )
    mtext("Read pairs assigned to cells (%)", side = 1, line = -2, cex = 1.5)

    # Draw box
    box(which = "outer", lty = "solid", col = "grey")
  })

  # Sequencing log-log plot
  output$seq_plot3 <- renderPlot({
    logLogPlot(sce()) +
      theme(
        title = element_text(size = 20),
        text = element_text(size = 16)
      )
  })

  # Panel plot 1
  output$panel_plot1 <- renderPlot({
    req(input$upload)
    sce_obj <- sce()
    metadata <- sce_obj@metadata

    plot(0, type = "n", axes = FALSE, ann = FALSE)
    mtext(metadata[["panel_name"]], side = 3, line = -2, cex = 3, 
          col = "#22A884FF")
    mtext("Panel used", side = 3, line = -4, cex = 1.5)

    # Print ean mapped reads per cell
    mtext(metadata[["n_amplicons"]], side = 1, line = -4, cex = 3, 
          col = "#22A884FF")
    mtext("Number of amplicons", side = 1, line = -2, cex = 1.5)
    box(which = "outer", lty = "solid", col = "grey")
  })

  # Panel plot 2
  output$panel_plot2 <- renderPlot({
    req(input$upload)
    sce_obj <- sce()
    metadata <- sce_obj@metadata

    genes <- vapply(str_split(rowData(sce)$id, "_"), function(x) x[3], 
                    character(1))
    plot(0, type = "n", axes = FALSE, ann = FALSE)
    mtext(length(unique(genes)), side = 3, line = -2, cex = 3, col = "#F66D7A")
    mtext("Number of Genes covered", side = 3, line = -4, cex = 1.5)

    # Print ean mapped reads per cell
    mtext(
      floor(as.numeric(metadata[["n_read_pairs_mapped_to_cells"]]) /
        as.numeric(metadata[["n_cells"]]) /
        as.numeric(metadata[["n_amplicons"]])),
      side = 1, line = -4,
      cex = 3, col = "#22A484FF"
    )
    mtext("Average read pairs per \namplicon and cell",
      side = 1, line = -0.5,
      cex = 1.5
    )
    box(which = "outer", lty = "solid", col = "grey")
  })

  # Panel karyoplot
  output$panel_plot3 <- renderPlot({
    plotAmpliconDistribution(sce = sce()) +
      theme(
        title = element_text(size = 20),
        text = element_text(size = 16)
      )
  })

  output$panel_plot4 <- renderPlot({
    plotNormalizedReadCounts(sce = sce()) +
      theme(
        title = element_text(size = 20),
        text = element_text(size = 16)
      )
  })

  output$panel_plot5 <- renderPlotly({
    plotPanelUniformity(sce = sce(), interactive = FALSE) +
      theme(
        title = element_text(size = 20),
        text = element_text(size = 16)
      )
  })

  # Tables ---------------------------------------------------------------------
  output$data_table_sample <- renderDataTable({
    sce_obj <- sce()
    metadata <- sce_obj@metadata
    metadata <- (metadata %>%
      unlist() %>%
      as.data.frame() %>%
      rownames_to_column())[c(16, 14, 6, 5, 4, 15), ] %>%
      as.data.frame() %>%
      mutate(rowname = gsub("n_", "Number of ", rowname)) %>%
      mutate(rowname = gsub("_", " ", rowname)) %>%
      mutate(rowname = str_to_title(rowname)) %>%
      datatable(.,
        extensions = "Buttons",
        options = list(
          pageLength = 6,
          dom = "Bt",
          buttons = list(
            list(extend = "csv", filename = 
                   paste0("scafari_sequencing_sample_", 
                          metadata[["sample_name"]])),
            list(extend = "excel", filename = 
                   paste0("scafari_sequencing_sample_", 
                          metadata[["sample_name"]])),
            list(extend = "pdf", filename = 
                   paste0("scafari_sequencing_sample_", 
                          metadata[["sample_name"]])),
            list(extend = "copy", filename = 
                   paste0("scafari_sequencing_sample_", 
                          metadata[["sample_name"]]))
          )
        ),
        colnames = rep("", ncol(.)), rownames = FALSE
      )
  })

  output$data_table_sequencing <- renderDataTable({
    sce_obj <- sce()
    metadata <- sce_obj@metadata
    rbind(
      "Total read pairs" = c(paste((as.numeric(metadata[["n_read_pairs"]])))),
      "Read pairs trimmed" = 
        c(paste((as.numeric(metadata[["n_read_pairs_trimmed"]])))),
      "Read pairs with valid barcodes" = 
        c(paste(round(as.numeric(metadata[["n_read_pairs_valid_cell_barcodes"]]))))
    ) %>%
      datatable(.,
        extensions = "Buttons",
        options = list(
          pageLength = 5, width = "100%",
          dom = "Bt",
          buttons = list(
            list(extend = "csv", filename = paste0("scafari_sequencing_overview_", metadata[["sample_name"]])),
            list(extend = "excel", filename = paste0("scafari_sequencing_overview_", metadata[["sample_name"]])),
            list(extend = "pdf", filename = paste0("scafari_sequencing_overview_", metadata[["sample_name"]])),
            list(extend = "copy", filename = paste0("scafari_sequencing_overview_", metadata[["sample_name"]]))
          )
        ), colnames = NULL
      )
  })

  output$data_table_mapping <- renderDataTable({
    sce_obj <- sce()
    metadata <- sce_obj@metadata
    rbind(
      "Reads mapped to genome (%)" = c(round((as.numeric(metadata[["n_reads_mapped"]]) / (as.numeric(metadata[["n_read_pairs"]]) * 2) * 100), digits = 2)),
      "Reads mapped to target (%)" = c(round((as.numeric(metadata[["n_reads_mapped_insert"]]) / (as.numeric(metadata[["n_read_pairs"]]) * 2) * 100), digits = 2))
    ) %>%
      datatable(.,
        extensions = "Buttons",
        options = list(
          pageLength = 5,
          width = "100%",
          dom = "Bt",
          buttons = list(
            list(extend = "csv", filename = paste0("scafari_sequencing_mapping", 
                                                   metadata[["sample_name"]])),
            list(extend = "excel", filename = paste0("scafari_sequencing_mapping", 
                                                     metadata[["sample_name"]])),
            list(extend = "pdf", filename = paste0("scafari_sequencing_mapping",
                                                   metadata[["sample_name"]])),
            list(extend = "copy", filename = paste0("scafari_sequencing_mapping", 
                                                    metadata[["sample_name"]]))
          )
        ),
        class = "display",
        colnames = rep("", ncol(.))
      )
  })

  ## Occurence of genes in panel -----------------------------------------------
  output$data_table_overview <- renderDataTable({
    annotateAmplicons(sce = sce())
  })
}
