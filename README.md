# scafari: Exploring scDNA-seq data

scafari is a Shiny app for the analysis of scDNA-seq data provided in `.h5` file. The analysis is divided into four steps, represented
by corresponding tabs: “Sequencing”, “Panel”, “Variants” and “Explore Variants”. 

# Requirements

To run scafari, you need R (Version 4.3.0 or higher) and R Shiny.


# Installation

To install scafari, download the repository. 

`git clone https://github.com/sophiewind/scafari.git`

# Running scafari
- Change directory to scafari: `cd scafari/`
- If you run scafari for the first time run setup.R to install the necessary packages: `Rscript setup.R`
- Now scafari is ready to run: `R -e "shiny::runApp('./app.R')"`
  - The libraries are loaded in the beginning of the run, this can take some time
  - When the libraries are loaded following line appears: `Listening on http://127.0.0.1:xxxx`. Click on the link, than your browser will connect to the shiny app

# Input

`.h5` files are accepted as input. After upload it is evaluated if all important information is inlcuded in the `.h5` file. 


# Features

## Upload
- data upload
- input checkup

## Sequencing
- basic information
- log-log plot
- sequencing statistics table
- mapping statistics table

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

## Panel
- hit the "Panel" tab


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
