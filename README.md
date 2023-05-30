<!-- badges: start -->
  [![Codecov test coverage](https://codecov.io/gh/dputhier/scigenex/branch/main/graph/badge.svg)](https://app.codecov.io/gh/dputhier/scigenex?branch=main)
  <!-- badges: end -->

# scigenex

WARNING: scigenex library is still in beta version.
WARNING: Online tutorial will be updated soon.

## Installation

### In the terminal

     R CMD INSTALL scigenex
     
     # In R
     
     library(scigenex)
     
### From R

The scigenex library is currently not available in CRAN or Bioc. To install from github, use:

     library(devtools)
     install_github("dputhier/scigenex")

To install from a tar in the terminal

    tar xvfz scigenex.tar.gz
    R CMD INSTALL dbfmcl
    
Or using the R interpreter:

    library(devtools)
    install("/path/to/the/package")
    library("scigenex")

## Example


		library(Seurat)
		library(scigenex)
		set_verbosity(1)

		# Load a dataset
		load_example_dataset("7871581/files/pbmc3k_medium")

		# Select informative genes
		res <- select_genes(pbmc3k_medium,
		                     distance = "pearson",
		                     row_sum=5)
		                     
		# Cluster informative features
		 
		## Construct and partition the graph
		res <- gene_clustering(res, 
		                       inflation = 1.5, 
		                       threads = 4)
		                        
		# Display the heatmap of gene clusters
		res <- top_genes(res)
		plot_heatmap(res, cell_clusters = Seurat::Idents(pbmc3k_medium))


### Documentation

Documentation (in progress) is available at [https://dputhier.github.io/scigenex/](https://dputhier.github.io/scigenex/).
