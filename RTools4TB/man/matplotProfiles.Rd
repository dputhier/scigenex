\name{matplotProfiles}
\alias{matplotProfiles}

\title{An internal plotting function based on the matplot function.}
\description{An internal function based on the \code{\link[graphics]{matplot}} function.}
\usage{
matplotProfiles(data, imgName, saveHTML = FALSE, X11 = TRUE, verbose = FALSE)
}

\arguments{
  \item{data}{a matrix}
  \item{imgName}{a prefix for image name.}
  \item{saveHTML}{a logical indicating if results should be stored in a HTML file.}
  \item{X11}{if set to TRUE a new graphic windows is generated.}
  \item{verbose}{if set to TRUE the function runs verbosely.}
}
\seealso{Used in the \code{\link{plotGeneExpProfiles}} function.}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\keyword{hplot}
