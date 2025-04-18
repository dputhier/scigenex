% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{cluster_set_from_seurat}
\alias{cluster_set_from_seurat}
\title{Transform a Seurat objects / FindAllMarkers result into a ClusterSet.}
\usage{
cluster_set_from_seurat(
  object = NULL,
  markers = NULL,
  which_slot = c("data", "counts", "sct"),
  assay_type = c("RNA", "Spatial")
)
}
\arguments{
\item{object}{A Seurat object.}

\item{markers}{A Seurat::FindAllMarkers() result or a named vector (clusters with gene_names as named).}

\item{which_slot}{One of 'data', 'counts' or 'sct'. The slot to extract from the seurat object to perform clustering analysis.
SCT is the recommended method from Seurat package when working with spatial transcriptomics data.}

\item{assay_type}{The type of assay ("RNA" or "Spatial").}
}
\description{
Transform a Seurat objects into a ClusterSet.
}
\examples{
## From a scRNA-seq/Seurat object
library(SeuratObject)
library(Seurat)
data("pbmc_small", package="SeuratObject")
cs <- cluster_set_from_seurat(pbmc_small, Seurat::FindAllMarkers(pbmc_small))
plot_heatmap(cs)
plot_heatmap(cs)
plot_heatmap(cs[1,])
plot_heatmap(cs, cell_clusters = Seurat::Idents(pbmc_small))
plot_heatmap(cs[1,Idents(pbmc_small) == "0"], 
             cell_clusters = Seurat::Idents(pbmc_small), label_size = 6)
plot_profiles(cs, ident = Seurat::Idents(pbmc_small))
}
