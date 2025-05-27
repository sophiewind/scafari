# scafari: Exploring scDNA-seq data

scafari is a Shiny app for the analysis of scDNA-seq data provided in `.h5` file. The analysis is divided into four steps, represented
by corresponding tabs: “Sequencing”, “Panel”, “Variants” and “Explore Variants”. 

# Requirements

To run scafari, you need R (Version 4.4 or higher) and R Shiny.

You can use bioconductor/bioconductor_docker (https://hub.docker.com/r/bioconductor/bioconductor_docker).

# Installation

To install scafari, open R and install the package.

```
install.packages("devtools")
devtools::install_github("sophiewind/scafari")
```

scafari is available as R package and as shiny app. A detailed description how to use the R packages is in the vignette in `/vignettes` directory.

# Running the scafari shiny app

- launch the app using `launchScafariShiny()`


# Input

- `.h5` files are accepted as input. After upload it is evaluated if all important information is inlcuded in the `.h5` file. 
- a test dataset is available in this repo in the `testdata` directory


# Features

## Upload
- data upload
- input checkup

## Sequencing
- basic information
- log-log plot
- sequencing statistics table
- mapping statistics table
- tapestri table
- R1 to R2 bar plot

## Panel
- panel information table
- amplicon overview table
- amplicon distribution plot
- amplicon performance plot
- amplicon uniformity plot

## Variants
- interactive filtering
- variant annotation
- filtered and annotated variants table
- VAF heatmap with genotype annotation
- genotype distribution and quality plot

## Explore variants
- interactive variant selection
- elbow plot (k-means)
- cluster plot (k-means)
- VAF, cluster and genotype plots splitted by clones

# Demo
- tables are downloaded to your `Downloads` directory.

## Upload
- just upload your `.h5` file

## Sequencing
- hit the "Sequencing" tab
- scroll through the page
- explore the interactive bar plot

## Panel
- hit the "Panel" tab
- the read counts are normalized and annotated in the background. Therefore, it 
can take some time until everything is processed and the plots are loaded 


## Variants
- click on the “Variants” tab 
- if necessary change default filtering parameters
- start with filtering the data by click on the `Apply Filtering`  button
- the filtering can take some time (☕)
- explore the filtered variants


## Explore variants
1. Variant selection
- hit the "Explore variants" tab
- with the results in the "Variants" tab you identify your variants of interest
- select those variants in the table on the left by clicking on them. You can deselect variants by clicking them a second time. On the right is a VAF heatmap which might help you to find you variants of interest
- if you selected all variants of interest click on `Select variants`
- the heatmap on the right will now be actualized and you can check you selection
  - if you want to change something you can change the selection on the table on the left and than hit `Select variants` again
  - if you want to continue hit `Continue with this variant selection`

2. Cell clustering
- select the number of centroids for the k-means clustering. The elbow plots will help you doing it
- if you know the number of clusteres you want to use. Change the number in `Number of clusters` by clicking on the small triangles
- hit `Start kmeans`

3. Explore clusters
- you can see which variants are included in the analysis under "Variants included in Clustering:"

# R Session
```
> sessionInfo()
R version 4.3.3 (2024-02-29)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 24.04.1 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0

locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8       
 [4] LC_COLLATE=de_DE.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8   
 [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                  LC_ADDRESS=C              
[10] LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
[1] grid      stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] markdown_1.13           factoextra_1.0.7        shinycustomloader_0.9.0 shinyBS_0.61.1         
 [5] shinyjs_2.1.0           org.Hs.eg.db_3.18.0     AnnotationDbi_1.64.1    Biobase_2.62.0         
 [9] clusterProfiler_4.10.1  biomaRt_2.58.2          reshape2_1.4.4          stringr_1.5.1          
[13] karyoploteR_1.28.0      regioneR_1.34.0         ComplexHeatmap_2.18.0   rhdf5_2.46.1           
[17] tibble_3.2.1            GenomicRanges_1.54.1    GenomeInfoDb_1.38.8     IRanges_2.36.0         
[21] S4Vectors_0.40.2        BiocGenerics_0.48.1     ggplot2_3.5.1           waiter_0.2.5           
[25] dplyr_1.1.4             DT_0.33                 shinycssloaders_1.1.0   shiny_1.9.1            

loaded via a namespace (and not attached):
  [1] splines_4.3.3               later_1.3.2                 BiocIO_1.12.0              
  [4] ggplotify_0.1.2             bitops_1.0-9                filelock_1.0.3             
  [7] polyclip_1.10-7             XML_3.99-0.17               rpart_4.1.23               
 [10] lifecycle_1.0.4             doParallel_1.0.17           MASS_7.3-60.0.1            
 [13] lattice_0.22-5              ensembldb_2.26.0            crosstalk_1.2.1            
 [16] backports_1.5.0             magrittr_2.0.3              sass_0.4.9                 
 [19] Hmisc_5.2-0                 rmarkdown_2.29              jquerylib_0.1.4            
 [22] yaml_2.3.10                 httpuv_1.6.15               cowplot_1.1.3              
 [25] DBI_1.2.3                   RColorBrewer_1.1-3          abind_1.4-8                
 [28] zlibbioc_1.48.2             purrr_1.0.2                 AnnotationFilter_1.26.0    
 [31] ggraph_2.2.1                biovizBase_1.50.0           RCurl_1.98-1.16            
 [34] yulab.utils_0.1.8           nnet_7.3-19                 tweenr_2.0.3               
 [37] VariantAnnotation_1.48.1    rappdirs_0.3.3              circlize_0.4.16            
 [40] GenomeInfoDbData_1.2.11     enrichplot_1.22.0           ggrepel_0.9.6              
 [43] tidytree_0.4.6              codetools_0.2-19            DelayedArray_0.28.0        
 [46] ggforce_0.4.2               DOSE_3.28.2                 xml2_1.3.6                 
 [49] tidyselect_1.2.1            shape_1.4.6.1               aplot_0.2.3                
 [52] farver_2.1.2                viridis_0.6.5               matrixStats_1.4.1          
 [55] BiocFileCache_2.10.2        base64enc_0.1-3             bamsignals_1.34.0          
 [58] GenomicAlignments_1.38.2    jsonlite_1.8.9              GetoptLong_1.0.5           
 [61] tidygraph_1.3.1             Formula_1.2-5               iterators_1.0.14           
 [64] foreach_1.5.2               tools_4.3.3                 progress_1.2.3             
 [67] treeio_1.26.0               Rcpp_1.0.13-1               glue_1.8.0                 
 [70] gridExtra_2.3               SparseArray_1.2.4           xfun_0.49                  
 [73] qvalue_2.34.0               MatrixGenerics_1.14.0       withr_3.0.2                
 [76] fastmap_1.2.0               rhdf5filters_1.14.1         fansi_1.0.6                
 [79] digest_0.6.37               gridGraphics_0.5-1          R6_2.5.1                   
 [82] mime_0.12                   colorspace_2.1-1            GO.db_3.18.0               
 [85] dichromat_2.0-0.1           RSQLite_2.3.7               utf8_1.2.4                 
 [88] tidyr_1.3.1                 generics_0.1.3              data.table_1.16.2          
 [91] rtracklayer_1.62.0          graphlayouts_1.2.0          prettyunits_1.2.0          
 [94] httr_1.4.7                  htmlwidgets_1.6.4           S4Arrays_1.2.1             
 [97] scatterpie_0.2.4            pkgconfig_2.0.3             gtable_0.3.6               
[100] blob_1.2.4                  XVector_0.42.0              shadowtext_0.1.4           
[103] htmltools_0.5.8.1           fgsea_1.28.0                ProtGenerics_1.34.0        
[106] clue_0.3-65                 scales_1.3.0                png_0.1-8                  
[109] ggfun_0.1.7                 knitr_1.48                  rstudioapi_0.17.1          
[112] rjson_0.2.23                nlme_3.1-164                checkmate_2.3.2            
[115] curl_6.0.0                  cachem_1.1.0                GlobalOptions_0.1.2        
[118] parallel_4.3.3              HDO.db_0.99.1               foreign_0.8-86             
[121] restfulr_0.0.15             pillar_1.9.0                vctrs_0.6.5                
[124] promises_1.3.0              dbplyr_2.5.0                xtable_1.8-4               
[127] cluster_2.1.6               htmlTable_2.4.3             evaluate_1.0.1             
[130] GenomicFeatures_1.54.4      cli_3.6.3                   compiler_4.3.3             
[133] bezier_1.1.2                Rsamtools_2.18.0            rlang_1.1.4                
[136] crayon_1.5.3                plyr_1.8.9                  fs_1.6.5                   
[139] stringi_1.8.4               viridisLite_0.4.2           BiocParallel_1.36.0        
[142] munsell_0.5.1               Biostrings_2.70.3           lazyeval_0.2.2             
[145] GOSemSim_2.28.1             Matrix_1.6-5                BSgenome_1.70.2            
[148] patchwork_1.3.0             hms_1.1.3                   bit64_4.5.2                
[151] Rhdf5lib_1.24.2             KEGGREST_1.42.0             SummarizedExperiment_1.32.0
[154] fontawesome_0.5.2           igraph_2.1.1                memoise_2.0.1              
[157] bslib_0.8.0                 ggtree_3.10.1               fastmatch_1.1-4            
[160] bit_4.5.0                   gson_0.1.0                  ape_5.8
```

Further tested on bioconductor/bioconductor_docker. BiocManager::version() >= '3.19'

