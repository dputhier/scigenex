# scigenex

WARNING: scigenex library is still in beta version. 

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

### Quick example on artificial data

      library(scigenex)
      set.seed(123)
      m <- matrix(rnorm(40000), nc=20)
      m[1:100,1:10] <- m[1:100,1:10] + 4
      m[101:200,11:20] <- m[101:200,11:20] + 3
      m[201:300,5:15] <- m[201:300,5:15] + -2
      res <- find_gene_clusters(data=m,
             distance_method="pearson",
             av_dot_prod_min = 0,
             inflation = 2,
             k=25,
             fdr = 10)
      plot_heatmap(res, ceil = 3, floor = -3)
      saveRDS(object = res, file = "ClusterSet_obj.RDS")
      
### Documentation

Documentation (in progress) is available at [https://dputhier.github.io/scigenex/](https://dputhier.github.io/scigenex/).
