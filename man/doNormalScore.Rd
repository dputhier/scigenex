\name{doNormalScore}
\alias{doNormalScore}

\title{Normal score transformation of a matrix.}
\description{
This function performs normal score transformation of a matrix. The normal score transformation ranks each column from lowest to the highest values and matches these ranks to equivalent ranks from a normal distribution.}
\usage{
doNormalScore(sdata,set.seed = 123)
}

\arguments{
  \item{sdata}{a matrix.}
  \item{set.seed}{specify seeds for random number generator.}
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
m[,1] <- m[,1] + 4
boxplot(as.data.frame(m))
m <- doNormalScore(m)
boxplot(as.data.frame(m))
}

\keyword{manip}
