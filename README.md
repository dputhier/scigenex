# dbfmcl

## Installation

### In the terminal

     R CMD INSTALL dbfmcl
     
     # In R
     
     library(dbfmcl)
     
### From R

To install from a private repo, generate a personal access token (PAT) in https://github.com/settings/tokens and supply to the *auth_token* argument. 

     library(devtools)
     install_github("JulieBvs/dbfmcl", auth_token="...", ref="develop")

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
               clustering=TRUE,
               k=25)
    plot_clust(res)
    write_clust(res, "ALL.sign.txt")
    
### With a subset of a seurat object

The Read10X()function reads in the output of the cellranger pipeline from 10X, returning a unique molecular identified (UMI) count matrix.
Data can be found here: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz 

    library(Seurat)
    library(dbfmcl)

    setwd("~/Downloads/")
    pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
    pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
    pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
    pbmc_df <- as.data.frame(pbmc@assays$RNA@data)
    dbf <- DBFMCL(pbmc_df, clustering = TRUE, k=30)
    p <- plot_clust(dbf)
    # Very long...
    # saveWidget(ggplotly(p), file = "test.html");
    # Faster 
    ggsave(p, height=max(dbf@cluster), filename="test.png")