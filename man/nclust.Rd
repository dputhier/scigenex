% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ClusterSet_class.R
\name{nclust}
\alias{nclust}
\title{Number of clusters in a ClusterSet object.
nclust}
\usage{
nclust(x)
}
\arguments{
\item{x}{The ClusterSet object}
}
\description{
The number of clusters in a ClusterSet object.
}
\examples{
load_example_dataset('7871581/files/pbmc3k_medium_clusters')
n_clust <- nclust(pbmc3k_medium_clusters)
}
\keyword{internal}
