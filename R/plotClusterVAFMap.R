plotClusterVAFMap <- function(sce, variants.of.interest, gg.clust){
  vaf.matrix.filtered <-  as.data.frame(t(assay(altExp(sce, 'variants'), 'VAF')))
  colnames(vaf.matrix.filtered) <- paste0(rowData(altExp(sce, 'variants'))$Gene, ':', rowData(altExp(sce, 'variants'))$id)
  rownames(vaf.matrix.filtered) <- paste0('cell', rownames(vaf.matrix.filtered)) # TODO cell ids
  vaf.matrix.filtered <- vaf.matrix.filtered[,variants.of.interest] 
  
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
  merged.var$variant <- factor(merged.var$variant, levels = sort(variants.of.interest))
  
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
  
  
  p <- ggplot(merged.var, aes(x = x, y = y, color = VAF)) +
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
      facet_grid(~factor(variant, levels = sort(variants.of.interest))) +
      theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        title = element_text(size = 20),
        text = element_text(size = 18)) +
      labs(x = gg.clust$labels$x, y = gg.clust$labels$y)
  return(p)
}
