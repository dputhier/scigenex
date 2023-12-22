<!-- README.md is generated from README.Rmd using devtools::build_readme(). Please edit that file -->

<img src="https://github.com/dputhier/scigenex/raw/master/inst/sticker/scigenex_logo.png" width="150" height="150" align="right"/>

    ## ✔ Setting active project to '/Users/puthier/Documents/git/project_dev/scigenex'

[![](https://img.shields.io/badge/lifecycle-stable-green.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Project Status: Active - The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![License:
MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://cran.r-project.org/web/licenses/MIT)
[![](https://img.shields.io/github/last-commit/dputhier/scigenex.svg)](https://github.com/dputhier/scigenex/commits/main)

# SciGeneX repository

## :arrow\_double\_down: Installation

### System requirements

The partitioning steps are currently performed using a system call to
the Markov Cluster (MCL) algorithm that presently limits the use of
DBF-MCL to unix-like platforms. Importantly, the `mcl` command should be
in your PATH and reachable from within R (see dedicated section).

### Step 1 - Installation of SciGeneX

#### From R

The scigenex library is currently not available in CRAN or Bioc. To
install from github, use:

    devtools::install_github("dputhier/scigenex")
    library(scigenex)

#### From the terminal

Download the *tar.gz* from github or clone the main branch. Uncompress
and run the following command from within the uncompressed scigenex
folder:

    R CMD INSTALL .

Then load the library from within R.

    library(scigenex)

### Step 2 - Installation of MCL

You may skip this step as the latest versions of SciGeneX will call
`scigenex::install_mcl()`to install MCL in `~/.scigenex` directory if
this program is not found in the PATH.

#### Installation of MCL using install\_mcl()

The `install_mcl()` has been developed to ease MCL installation. This
function should be call automatically from within R when calling the
`gene_clustering()` function. If `install_mcl()` does not detect MCL in
the PATH it will install it in `~/.scigenex`.

#### Installation of MCL from source

One also can install MCL from source using the following code.

    # Download the latest version of mcl 
    wget http://micans.org/mcl/src/mcl-latest.tar.gz
    # Uncompress and install mcl
    tar xvfz mcl-latest.tar.gz
    cd mcl-xx-xxx
    ./configure
    make
    sudo make install
    # You should get mcl in your path
    mcl -h

#### Installation of MCL from sources

Finally you may install MCL using conda. Importantly, the mcl command
should be available in your PATH from within R.

    conda install -c bioconda mcl

## Example

The scigenex library contains several datasets including the
pbmc3k\_medium which is a subset from pbmc3k 10X dataset.

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

## :book: Documentation

Documentation (in progress) is available at
<https://dputhier.github.io/scigenex/>.
