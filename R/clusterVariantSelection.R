#' Function: clusterVariantSelection
#' This function takes selected variants and performs clustering on them.
#'
#' @param sce A SingleCellExperiment object containing the single-cell data on
#' which clustering will be performed.
#' @param variants.of.interest A vector or list specifying the variants of
#' interest to be selected for clustering.
#' @param n.clust An integer specifying the number of clusters.
#' @param resolution The resolution parameter to use. Higher resolutions lead
#' to more smaller communities, while lower resolutions lead to fewer larger
#' communities.
#' @param eps.value Size (radius) of the epsilon neighborhood. Can be omitted
#' if x is a frNN object.
#' @param min.pts Number of minimum points required in the eps neighborhood
#' for core points (including the point itself). By default cell number/100.
#' @param method Clustering method. Either k-means, dbscan or leiden.
#' @return A list with clustering results and a ggplot-object.
#'
#' @references https://cran.r-project.org/web/packages/dbscan/readme/
#' README.html#ref-hahsler2019dbscan
#'
#' @examples
#' # Assume `sce` is a SingleCellExperiment object with variants in altExp()
#' sce_filtered <- readRDS(system.file("extdata", "sce_filtered_demo.rds",
#'     package = "scafari"
#' ))
#' clusterplot <- clusterVariantSelection(
#'     sce = sce_filtered,
#'     variants.of.interest = c(
#'         "FLT3:chr13:28610183:A/G",
#'         "KIT:chr4:55599436:T/C",
#'         "TP53:chr17:7577427:G/A",
#'         "TET2:chr4:106158216:G/A"
#'     ),
#'     n.clust = 4
#' )
#'
#' @export
clusterVariantSelection <- function(sce, variants.of.interest, n.clust,
                                    method = "k-means",
                                    eps.value = 0.2,
                                    resolution = NULL,
                                    min.pts = NULL) {
    # Check that the input is a SingleCellExperiment object
    if (!inherits(sce, "SingleCellExperiment")) {
        stop("The input must be a SingleCellExperiment object.")
    }

    # Check that 'variants' altExp exists and contains a 'VAF' assay
    if (!"variants" %in% altExpNames(sce)) {
        stop(
            "The SingleCellExperiment object must contain 'variants' as an ",
            "alternate experiment."
        )
    }
    if (!"VAF" %in% assayNames(altExp(sce, "variants"))) {
        stop("The 'variants' alternate experiment must contain a 'VAF' assay.")
    }

    # Check that variants.of.interest is non-empty and exists in the data
    if (length(variants.of.interest) == 0) {
        stop("variants.of.interest must be a non-empty vector.")
    }

    vaf.matrix.filtered <- as.data.frame(t(assay(
        altExp(sce, "variants"),
        "VAF"
    )))
    colnames(vaf.matrix.filtered) <-
        paste0(
            rowData(altExp(sce, "variants"))$Gene, ":",
            rowData(altExp(sce, "variants"))$id
        )

    if (!all(variants.of.interest %in% colnames(vaf.matrix.filtered))) {
        stop("All variants.of.interest must exist in the VAF matrix columns.")
    }

    if (!method %in% c("leiden", "k-means", "dbscan")) {
        stop(
            "Invalid clustering method. Only 'leiden' and 'k-means' are ",
            "valid for this parameter."
        )
    }


    df <- vaf.matrix.filtered[, variants.of.interest] # selected_variants()]
    df <- na.omit(df)
    df <- scale(df)

    # k-means ------------------------------------------------------------------
    if (method == "k-means") {
        if (is.null(n.clust)) {
            stop(
                "Parameter 'n.clust' must be provided when using 'k-means' ",
                "clustering."
            )
        }

        # Determining Optimal Clusters
        kmeans_result <- kmeans(df, centers = n.clust, nstart = 25)

        # Generate the cluster plot and store it in shared_data
        gg.clust <- fviz_cluster(kmeans_result,
            data = df,
            ellipse.type = "norm",
            geom = "point"
        ) +
            theme_default()

        return(list(k_means = kmeans_result, clusterplot = gg.clust))

        # Leiden ---------------------------------------------------------------
    } else if (method == 'leiden'){
        tryCatch({
            # Perform PCA
            pca_result <- prcomp(df, center = TRUE, scale. = TRUE)
            pc_scores <- pca_result$x
            
            # Utilize RANN for KNN graph construction
            neighborhood <- RANN::nn2(pc_scores)
            
            # Create adjacency matrix for KNN graph
            adjacency_matrix <- matrix(0, ncol = nrow(pc_scores), 
                                    nrow = nrow(pc_scores))
            for (i in seq_len(nrow(pc_scores))) {
                adjacency_matrix[i, neighborhood$nn.idx[i, ]] <- 1
            }
            
            # Create graph and apply Leiden clustering
            knn_graph <- igraph::graph_from_adjacency_matrix(
                adjacency_matrix, mode = "undirected")
            leiden_results <- igraph::cluster_leiden(knn_graph, 
                                                    resolution = resolution)
            
            cluster <- leiden_results$membership
            
            # Handle scenario with too many clusters
            if (length(unique(cluster)) > 20) {
                return(list(
                    leiden = paste0("Error: Too many clusters found (",
                                    length(unique(cluster)), "). This value ",
                                    "is > 21. Adjust resolution."),
                    clusterplot = "Error: Too many clusters. Adjust resolution."
                ))
            }
            
            # Prepare cluster data frame for plotting
            cluster_data <- as.data.frame(cbind(pc_scores, cluster))
            cluster_data$cluster <- as.factor(cluster_data$cluster)
            colnames(cluster_data)[c(1, 2)] <- c('x', 'y')
            
            # Generate ggplot
            clust_plot <- ggplot(cluster_data, aes(x = x, y = y, 
                                                color = cluster)) +
                geom_point() +
                labs(
                    x = paste0('PC1 (', round(
                        summary(pca_result)$importance[2,1], 3) * 100, '%)'),
                    y = paste0('PC2 (', round(
                        summary(pca_result)$importance[2,2], 3) * 100, '%)')
                ) +
                theme_default()
            
            return(list(leiden_results = leiden_results, 
                        clusterplot = clust_plot))
            
        }, error = function(e) {
            # Catch any unexpected errors and return as structured messages
            return(list(
                leiden = "Unexpected error occurred during leiden clustering.",
                clusterplot = paste("Error: ", e$message)
            ))
        })

        # DBSCAN ---------------------------------------------------------------
    } else if (method == "dbscan") {
        # Set the parameters for DBSCAN
        if (is.null("min.pts")) min.pts <- dim(df)[1] / 100

        dbscan::kNNdistplot(df, k = min.pts)
        
        # Run DBSCAN
        dbscan_result <- dbscan(df, eps = eps.value, minPts = min.pts)
        
        pca_result <- prcomp(df, center = TRUE, scale. = TRUE)
        pca_scores <- pca_result$x
        cluster_data <- as.data.frame(cbind(pca_scores, dbscan_result$cluster))
        
        # Bring cluster from 0 to n to 1 to n + 1
        cluster_data[, ncol(cluster_data)] <- 
            cluster_data[, ncol(cluster_data)] + 1
        
        # Change colnames -> comparable between clustering methods
        colnames(cluster_data)[c(1, 2)] <- c('x', 'y')
        colnames(cluster_data)[ncol(cluster_data)] <- 'cluster'
        
        cluster_data$cluster <- as.factor(cluster_data$cluster)
        clust_plot <- ggplot(cluster_data, aes(x = x, y = y, color = cluster)) +
            geom_point() + 
            labs(x = paste0('PC1 (', round(summary(pca_result)$importance[2,1],
                                        3)*100, '%)'),
                y = paste0('PC1 (', round(summary(pca_result)$importance[2,2], 
                                        3)*100, '%)')) +
            stat_ellipse(aes(fill = cluster), 
                        alpha = .25, geom = 'polygon') +
            theme_default()
        return(list(dbscan = dbscan, clusterplot = clust_plot))
    }
}
