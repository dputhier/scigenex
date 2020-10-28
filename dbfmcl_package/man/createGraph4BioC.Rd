\name{createGraph4BioC}
\alias{createGraph4BioC}
\title{Creates a graph based on a request to the TrancriptomeBrowser database.}
\description{
Results from a request to TBrowserDB can be displayed as a graph using the createGraph4BioC function. Given a request (e.g., "XBP1 & ESR1 & GATA3"), createGraph4BioC retrieves the list of signature IDs that verify the constrain. A list of gene falling in at least one of the clusters is next computed. A matrix is created that will record for each gene the number of time they have been observed in the same signature (only genes falling in a significant proportion of signatures are conserved). This adjacency matrix is used to create a graph. 
}
\usage{
createGraph4BioC(request = NULL, prop = 50)
}

\arguments{
  \item{request}{a Boolean query with HUGO IDs (e.g., \code{"CD4 & CD3E & CD3D"})}
  \item{prop}{only genes falling in a significant proportion of signatures are conserved. For instance, \code{prop=50} indicates that only genes falling in at least 50 \% of the signatures will be conserved.}
}

\details{
In order to create an adjacence matrix, this function use the both \code{\link{getSignatures}} and \code{\link{getExpressionMatrix}} functions.

The "value" argument to this function may contain Boolean operators (see help section on TBrowser web site for more informations, \url{http://tagc.univ-mrs.fr/tbrowser} or the details section of \code{\link{getSignatures}}). 

\code{&}    : AND
\code{|}    : OR
\code{!}    : NOT , (used in conjonction with &)

}
\value{
an adjacence matrix.
}
\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\examples{
\dontrun{
# Create a graph based on expression signatures that contain "ESR1 & GATA3 & XBP1".
# Only genes observed in 80% of this signatures will be used.  
library(biocGraph)
adjMat <- createGraph4BioC(request="ESR1 & GATA3 & XBP1", prop=80)
g1 <- new("graphAM", adjMat=adjMat)
plot(g1,"fdp")
g2 <- as(g1,"graphNEL")
}
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\keyword{aplot}
\keyword{hplot}