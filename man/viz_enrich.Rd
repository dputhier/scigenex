% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/enrichment_analysis.R
\name{viz_enrich}
\alias{viz_enrich}
\title{Visualization of enrichment analyses results from a ClusterSet object}
\usage{
viz_enrich(
  object,
  clusters = "all",
  type = "dotplot",
  nb_terms = 20,
  terms_size = 50,
  font_size = 4
)
}
\arguments{
\item{object}{A \code{ClusterSet} object.}

\item{clusters}{A vector of cluster id to plot.}

\item{type}{A vector providing the type of plot.}

\item{nb_terms}{An integer indicating the number of terms in the plot.}

\item{terms_size}{An interger indicating the wrap length of terms}

\item{font_size}{The font size for the x axis.}
}
\value{
A \code{ClusterSet} object
}
\description{
Display enrichment analyses results for a ClusterSet object.
}
\keyword{internal}
