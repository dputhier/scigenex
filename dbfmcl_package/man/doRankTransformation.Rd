\name{doRankTransformation}
\alias{doRankTransformation}

\title{Rank transformation of a matrix.}
\description{
This function performs rank transformation of a matrix (each column is transformed so that each values is replaced by its corresponding ranks.).
}
\usage{
doRankTransformation(data = NULL)
}

\arguments{
  \item{data}{a matrix}
}

\value{
A matrix.
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\examples{
m <- matrix(rnorm(1000),nc=4)
m <- doRankTransformation(m)
head(m)
}

\keyword{manip}
