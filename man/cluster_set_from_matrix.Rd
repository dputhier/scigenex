% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert.R
\name{cluster_set_from_matrix}
\alias{cluster_set_from_matrix}
\title{Transform any matrix and list into a ClusterSet object.}
\usage{
cluster_set_from_matrix(object = NULL, markers = NULL)
}
\arguments{
\item{object}{A matrix or data.frame.}

\item{markers}{A list of vector containing the gene sets}
}
\description{
Transform any matrix (e.g expression matrix) and list (e.g markers obtained through kmeans 
or any partitioning algorithm)  into a ClusterSet object.
}
\examples{
m <- create_3_rnd_clust()[1:300,] 
rownames(m) <- paste0("gene", 1:300)
markers <- list(a=paste0("gene", 1:100), 
                b=paste0("gene", 101:200),
                c=paste0("gene", 201:300))
cs <- cluster_set_from_matrix(m, markers)
plot_heatmap(cs, interactive = FALSE)

}
