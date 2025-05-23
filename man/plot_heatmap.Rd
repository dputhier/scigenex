% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization.R
\name{plot_heatmap}
\alias{plot_heatmap}
\title{plot_heatmap}
\usage{
plot_heatmap(
  object,
  center = TRUE,
  ceil = 1,
  floor = -1,
  cell_clusters = NULL,
  show_dendro = TRUE,
  use_top_genes = FALSE,
  interactive = TRUE,
  name = NULL,
  xlab = NULL,
  ylab = NULL,
  colorbar_name = "Exp. level",
  show_legend = TRUE,
  colors = colors_for_gradient("Ju1"),
  colors_cell_clusters = NULL,
  row_labels = TRUE,
  col_labels = FALSE,
  label_size = 10,
  line_size_vertical = 3,
  line_size_horizontal = 3,
  link = c("average", "complete", "single", "ward.D", "ward.D2", "mcquitty")
)
}
\arguments{
\item{object}{A ClusterSet object.}

\item{center}{A logical to indicate whether to center row.}

\item{ceil}{A value for ceiling (NULL for no ceiling). Ceiling is performed centering.}

\item{floor}{A value for flooring (NULL for no flooring). Flooring is performed after centering.}

\item{cell_clusters}{A vector of cell clusters with cell barcodes as names.}

\item{show_dendro}{A logical to indicate whether to show column dendrogram.}

\item{use_top_genes}{A logical to indicate whether to use highly similar genes in the slot top_genes of ClusterSet.}

\item{interactive}{A logical to indicate if the heatmap should be interactive.}

\item{name}{A title for the heatmap (if interactive is TRUE).}

\item{xlab}{A title for the x axis (if interactive is TRUE).}

\item{ylab}{A title for the y axis (if interactive is TRUE).}

\item{colorbar_name}{A title for the colorbar.}

\item{show_legend}{A logical to indicate whether to show colorbar.}

\item{colors}{A vector of colors.}

\item{colors_cell_clusters}{A named vector of colors for cell identity annotations.}

\item{row_labels}{A logical to indicate whether to show row labels.}

\item{col_labels}{A logical to indicate whether to show col labels.}

\item{label_size}{A value for label font size.}

\item{line_size_vertical}{An integer for the size of horizontal white line which separate gene clusters.}

\item{line_size_horizontal}{An integer for the size of vertical white line  which separate cell clusters.}

\item{link}{The aggloremative criterion for hierarchical clustering. One of "average", "complete" or "single". 
Default to average.}
}
\value{
Iheatmap-class object.
}
\description{
Create a heatmap from a ClusterSet object.
}
\examples{
library(Seurat)
# Set verbosity to 1 to display info messages only.
set_verbosity(1)

# Load datasets
load_example_dataset('7871581/files/pbmc3k_medium_clusters')
load_example_dataset('7871581/files/pbmc3k_medium')

# Plot heatmap of 'top genes' of all gene clusters
pbmc3k_medium_clusters <- top_genes(pbmc3k_medium_clusters)
plot_heatmap(pbmc3k_medium_clusters, use_top_genes=TRUE)

# Plot heatmap of gene cluster 1
plot_heatmap(pbmc3k_medium_clusters[1,]) 

# Plot heatmap of gene cluster 1 and 3
plot_heatmap(pbmc3k_medium_clusters[c(1,3),])  

# Plot heatmap of 'top genes' of all gene clusters
# with cell ordered according to Seurat results
plot_heatmap(pbmc3k_medium_clusters, use_top_genes=TRUE, cell_clusters=Seurat::Idents(pbmc3k_medium))

# Plot heatmap of 'top genes' of all gene clusters
# with cell ordered according to Seurat results
# (non interactive version)
plot_heatmap(pbmc3k_medium_clusters, use_top_genes=TRUE, 
             cell_clusters=Seurat::Idents(pbmc3k_medium),
             interactive=FALSE, label_size = 2)
}
