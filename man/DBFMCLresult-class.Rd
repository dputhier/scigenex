\name{DBFMCLresult}
\docType{class}
\alias{DBFMCLresult}
\alias{DBFMCLresult}
\alias{show,DBFMCLresult-method}
\title{Class to store DBFMCL results.}
\description{This class represents the results of the \code{\link{DBFMCL}} algorithm.}
\section{Objects from the Class}{
Objects can be created by calls of the form
\code{  new('DBFMCLresult',}\cr
\code{    name          = ...., # Object of class character}\cr
\code{    data          = ...., # Object of class matrix}\cr
\code{    cluster       = ...., # Object of class vector}\cr
\code{    size          = ...., # Object of class vector}\cr
\code{    center        = ...., # Object of class matrix}\cr
\code{    parameters    = ...., # Object of class list}\cr
\code{  )}
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{Object of class \code{"character"}, an analysis identifiant (by default "exprs").}
    \item{\code{data}:}{Object of class \code{"matrix"}, a subset of the original matrix containing the coordinates of the selected elements.}
    \item{\code{cluster}:}{Object of class \code{"vector"}, a vector of integers indicating the cluster to which each point is allocated.}
    \item{\code{size}:}{Object of class \code{"vector"}, the number of points in each cluster.}
    \item{\code{center}:}{Object of class \code{"matrix"}, a matrix of cluster centres.}
    \item{\code{parameters}:}{Object of class \code{"list"}, a list of all used DBFMCL parameters: normalizationMethod, distanceMethod, k, random, fdr, set.seed, inflation.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "DBFMCLresult")}, provides informations about the object.}
    \item{plot_dbf}{\code{signature(object = "DBFMCLresult")}, plot a DBFMCLresult.}
    \item{write_dbf}{\code{signature(object = "DBFMCLresult")}, write a DBFMCLresult.}
	 }
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}
\author{Bergon A., Lopez F., Textoris J. and Puthier D.}



\seealso{\code{\link{DBFMCL}}}
\examples{
obj <- new("DBFMCLresult")
obj
}
\keyword{classes}
