\name{DBF}
\alias{DBF}

\title{Density-Based Filtering.}
\description{
This function is an internal function used by \code{\link{DBFMCL}} to detect informative elements (\emph{i.e.}, those that belong to dense regions). User should not use this function. Instead they can use the \code{\link{DBFMCL}} function with \code{clustering} argument set to \code{FALSE}.}

\usage{
DBF(data, name = NULL, distance.method = c("spearman", "pearson", "euclidean"), silent = FALSE, k = 100, random = 3, fdr = 10, memory.used = 1024, set.seed = 123, returnRank = FALSE)
}

\arguments{
  \item{data}{a matrix or data.frame}
  \item{name}{a prefix for the file name}
 \item{distance.method}{a method to compute the distance to the k-th nearest neighbor. One of "pearson" (Pearson's correlation coefficient-based distance), "spearman" (Spearman's rho-based distance) or "euclidean".}
  \item{silent}{if set to TRUE (default), the progression of distance matrix calculation is not displayed.}
  \item{k}{the neighborhood size.}
  \item{random}{the number of simulated distributions S to compute. By default \code{random = FALSE}.}
  \item{fdr}{a value for the false discovery rate.}
  \item{memory.used}{size of the memory used to store part of the distance matrix. The subsequent sub-matrix is used to computed simulated distances to the k-th nearest neighbor (see detail section).}
  \item{set.seed}{specify seeds for random number generator.}  
  \item{returnRank}{This argument modifies the output. Given a set of elements conserved after the filtering step of the DBFMCL algorithm, if \code{returnRank = TRUE} their expression values are replaced by their corresponding ranks in the input matrix.}
}
\details{
See \code{\link{DBFMCL}}
}

\section{Warnings}{Works only on UNIX-alikes platforms.}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}
\seealso{\code{\link{DBFMCL}}, \code{\link{createSignatures4TB}}}
\keyword{manip}
