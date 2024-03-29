% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_analysis.R
\name{enrich_go,ClusterSet-method}
\alias{enrich_go,ClusterSet-method}
\title{Enrichment analysis using GO database of gene clusters from a ClusterSet object}
\usage{
\S4method{enrich_go}{ClusterSet}(object, species = "Hsapiens", ontology = "all")
}
\arguments{
\item{object}{an object of class \code{ClusterSet}.}

\item{species}{a character string indicating the species name,
as a concatenation of the first letter of the name (uppercase) and the family name (lowercase),
 e.g human -> Hsapiens, mouse -> Mmusculus}

\item{ontology}{a character string indicating the ontology to use for the enrichment analysis. One of "BP", "MF", and "CC" ontology, or "ALL".}
}
\value{
A \code{ClusterSet} object
}
\description{
This function performs enrichment analysis of all gene clusters from a ClusterSet object
 and store the results in the cluster_annotations slot.
}
\examples{
load_example_dataset("7871581/files/pbmc3k_medium_clusters")
load_example_dataset("8028126/files/pbmc3k_medium_clusters_enr")
# pbmc3k_medium_clusters_enr <- enrich_go(pbmc3k_medium_clusters)
pbmc3k_medium_clusters_enr
is(pbmc3k_medium_clusters_enr@gene_cluster_annotations)
}
