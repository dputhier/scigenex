% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clusterset_report.R
\name{cluster_set_report}
\alias{cluster_set_report}
\title{Create a report from a ClusterSet and Seurat object. DEPRECATED}
\usage{
cluster_set_report(
  clusterset_object = NULL,
  seurat_object = NULL,
  file_path = NULL,
  force = FALSE,
  report_title = "Scigenex_report",
  report_author = "Unknown",
  heatmap_colors = colors_for_gradient("Ju1"),
  heatmap_color_ident = NULL,
  coord_flip = TRUE,
  pt_size = 2.75,
  go_info = c("ID", "Description", "GeneRatio", "pvalue", "qvalue", "Count"),
  verbosity = 1,
  pandoc_heap_size = "512m",
  report_id = rlang::hash(clusterset_object)
)
}
\arguments{
\item{clusterset_object}{The ClusterSet object.}

\item{seurat_object}{The Seurat object.}

\item{file_path}{A file path where to store the report (html extension).}

\item{force}{Whether to force erase output file.}

\item{report_title}{A title for the report.}

\item{report_author}{Names of report authors.}

\item{heatmap_colors}{Colors for the heatmap gradient.}

\item{heatmap_color_ident}{Colors for cell classes.}

\item{coord_flip}{Whether to flip spatial coordinates.}

\item{pt_size}{Size of the spots.}

\item{go_info}{Info about GO stats. Should be from the list: "ONTOLOGY", "ID", "Description", 
"GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count"}

\item{verbosity}{The verbosity level (0-3, defaut to 1.)}

\item{pandoc_heap_size}{Fix the maximum memory used by pandoc. See https://pandoc.org/MANUAL.html.}

\item{report_id}{Used to label code chunks.}
}
\description{
Create a report from a ClusterSet and Seurat object. DEPRECATED. Use scigenex_report() please.
}
\examples{
# Load datasets
set_verbosity(3)
library(Seurat)
load_example_dataset("7870305/files/lymph_node_tiny_clusters_2")
load_example_dataset("7870305/files/lymph_node_tiny_2")
# Create a report
cluster_set_report(clusterset_object = lymph_node_tiny_clusters_2[1,],
                   seurat_object=lymph_node_tiny_2,
                   pt_size=6)
}
