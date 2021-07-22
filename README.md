# dbfmcl

WARNING: dbfmcl library is still in beta version. 

## Installation

### In the terminal

     R CMD INSTALL dbfmcl
     
     # In R
     
     library(dbfmcl)
     
### From R

The dbfmcl library is currently not available in CRAN or Bioc. To install from github, use:

     library(devtools)
     install_github("dputhier/dbfmcl")

To install from a tar in the terminal

    tar xvfz dbfmcl.tar.gz
    R CMD INSTALL dbfmcl
    
Or using the R interpreter:

    library(devtools)
    install("/path/to/the/package")
    library("dbfmcl")

## Example

### Quick example on artificial data

      m <- matrix(rnorm(80000), nc=20)
      m[1:100,1:10] <- m[1:100,1:10] + 4
      m[101:200,11:20] <- m[101:200,11:20] + 3
      m[201:300,5:15] <- m[201:300,5:15] + -2
      res <- DBFMCL(data=m,
                    distance_method="pearson",
                    av_dot_prod_min = 0,
                    inflation = 1.2,
                    k=25,
                    fdr = 10)
      plot_clust(res, ceil = 10, floor = -10)
      write_clust(res, "ALL.sign.txt")
    
### Documentation

Documentation (in progress) is available at [https://dputhier.github.io/dbfmcl/](https://dputhier.github.io/dbfmcl/).