# scafari: Exploring scDNA-seq data

scafari is a Shiny app for the analysis of scDNA-seq data provided in `.h5` file. The analysis is divided into four steps, represented
by corresponding tabs: “Sequencing”, “Panel”, “Variants” and “Explore Variants”. 

# Requirements

To run scafari, you need R (Version 4.5 or higher) and R Shiny.

scafari was tested in the latest bioconductor devel docker (docker pull bioconductor/bioconductor_docker:devel).


# Installation

To install scafari, open R and install the package.

```
BiocManager::install('scafari')
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
R version 4.5.0 (2025-04-11)
Platform: x86_64-pc-linux-gnu
Running under: Ubuntu 24.04.2 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.12.0 
LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.12.0  LAPACK version 3.12.0

locale:
 [1] LC_CTYPE=de_DE.UTF-8       LC_NUMERIC=C               LC_TIME=de_DE.UTF-8        LC_COLLATE=de_DE.UTF-8     LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=de_DE.UTF-8    LC_PAPER=de_DE.UTF-8      
 [8] LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       

time zone: Europe/Berlin
tzcode source: system (glibc)

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] scafari_0.99.12             SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0         IRanges_2.42.0             
 [8] S4Vectors_0.46.0            BiocGenerics_0.54.0         generics_0.1.4              MatrixGenerics_1.20.0       matrixStats_1.5.0           rhdf5_2.52.0                shiny_1.10.0               

loaded via a namespace (and not attached):
  [1] fs_1.6.6                 ProtGenerics_1.40.0      bitops_1.0-9             devtools_2.4.5           fontawesome_0.5.3        httr_1.4.7               RColorBrewer_1.1-3       doParallel_1.0.17       
  [9] profvis_0.4.0            tools_4.5.0              backports_1.5.0          R6_2.6.1                 DT_0.33                  lazyeval_0.2.2           rhdf5filters_1.20.0      GetoptLong_1.0.5        
 [17] urlchecker_1.0.1         withr_3.0.2              prettyunits_1.2.0        GGally_2.2.1             gridExtra_2.3            textshaping_1.0.1        factoextra_1.0.7         cli_3.6.5               
 [25] shinyjs_2.1.0            ggbio_1.56.0             labeling_0.4.3           sass_0.4.10              askpass_1.2.1            systemfonts_1.2.3        Rsamtools_2.24.0         txdbmaker_1.4.1         
 [33] foreign_0.8-86           dbscan_1.2.2             R.utils_2.13.0           stringdist_0.9.15        dichromat_2.0-0.1        sessioninfo_1.2.3        styler_1.10.3            BSgenome_1.76.0         
 [41] rstudioapi_0.17.1        RSQLite_2.3.11           shape_1.4.6.1            BiocIO_1.18.0            crosstalk_1.2.1          car_3.1-3                dplyr_1.1.4              Matrix_1.6-5            
 [49] waldo_0.6.1              abind_1.4-8              R.methodsS3_1.8.2        lifecycle_1.0.4          yaml_2.3.10              carData_3.0-5            biocViews_1.76.0         SparseArray_1.8.0       
 [57] BiocFileCache_2.16.0     grid_4.5.0               blob_1.2.4               promises_1.3.2           crayon_1.5.3             miniUI_0.1.2             lattice_0.22-5           GenomicFeatures_1.60.0  
 [65] KEGGREST_1.48.0          sys_3.4.3                pillar_1.11.0            knitr_1.50               ComplexHeatmap_2.24.0    rjson_0.2.23             codetools_0.2-19         glue_1.8.0              
 [73] data.table_1.17.2        remotes_2.5.0            vctrs_0.6.5              png_0.1-8                testthat_3.2.3           gtable_0.3.6             cachem_1.1.0             xfun_0.52               
 [81] S4Arrays_1.8.0           mime_0.13                iterators_1.0.14         ellipsis_0.3.2           usethis_3.1.0            bit64_4.6.0-1            progress_1.2.3           filelock_1.0.3          
 [89] rprojroot_2.0.4          R.cache_0.17.0           bslib_0.9.0              rpart_4.1.23             colorspace_2.1-1         DBI_1.2.3                Hmisc_5.2-3              shinycustomloader_0.9.0 
 [97] nnet_7.3-19              tidyselect_1.2.1         processx_3.8.6           waiter_0.2.5             bit_4.6.0                compiler_4.5.0           curl_6.2.2               httr2_1.1.2             
[105] graph_1.86.0             BiocCheck_1.44.2         htmlTable_2.4.3          xml2_1.3.8               plotly_4.10.4            desc_1.4.3               DelayedArray_0.34.1      rtracklayer_1.68.0      
[113] checkmate_2.3.2          scales_1.4.0             RBGL_1.84.0              callr_3.7.6              rappdirs_0.3.3           stringr_1.5.1            digest_0.6.37            shinyBS_0.61.1          
[121] rmarkdown_2.29           XVector_0.48.0           htmltools_0.5.8.1        pkgconfig_2.0.3          base64enc_0.1-3          dbplyr_2.5.0             fastmap_1.2.0            ensembldb_2.32.0        
[129] rlang_1.1.6              GlobalOptions_0.1.2      htmlwidgets_1.6.4        UCSC.utils_1.4.0         farver_2.1.2             jquerylib_0.1.4          jsonlite_2.0.0           BiocParallel_1.42.0     
[137] R.oo_1.27.1              VariantAnnotation_1.54.1 RCurl_1.98-1.17          magrittr_2.0.3           Formula_1.2-5            GenomeInfoDbData_1.2.14  credentials_2.0.2        Rhdf5lib_1.30.0         
[145] Rcpp_1.0.14              shinycssloaders_1.1.0    stringi_1.8.7            brio_1.1.5               MASS_7.3-60.0.1          org.Hs.eg.db_3.21.0      plyr_1.8.9               pkgbuild_1.4.7          
[153] ggstats_0.9.0            ggrepel_0.9.6            parallel_4.5.0           Biostrings_2.76.0        hms_1.1.3                circlize_0.4.16          ps_1.9.1                 igraph_2.1.4            
[161] ggpubr_0.6.0             RUnit_0.4.33.1           markdown_2.0             ggsignif_0.6.4           reshape2_1.4.4           biomaRt_2.64.0           pkgload_1.4.0            XML_3.99-0.18           
[169] evaluate_1.0.3           biovizBase_1.56.0        BiocManager_1.30.25      foreach_1.5.2            httpuv_1.6.16            RANN_2.6.2               tidyr_1.3.1              openssl_2.3.2           
[177] purrr_1.0.4              clue_0.3-66              ggplot2_3.5.2            BiocBaseUtils_1.10.0     broom_1.0.8              xtable_1.8-4             restfulr_0.0.15          AnnotationFilter_1.32.0 
[185] roxygen2_7.3.2           rstatix_0.7.2            later_1.4.2              ragg_1.4.0               viridisLite_0.4.2        gert_2.1.5               OrganismDbi_1.50.0       tibble_3.2.1            
[193] memoise_2.0.1            AnnotationDbi_1.70.0     GenomicAlignments_1.44.0 cluster_2.1.6            BiocStyle_2.36.0            gson_0.1.0                  ape_5.8
```


