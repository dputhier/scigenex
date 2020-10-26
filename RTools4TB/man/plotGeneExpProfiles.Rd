\name{plotGeneExpProfiles}
\alias{plotGeneExpProfiles}

\title{Visualization of transcriptional signature profiles.}
\description{This function is used to visualize expression profiles of transcriptional signatures stored in a DBFMCLresult object or in a tab-delimited file.}

\usage{
plotGeneExpProfiles(data = NULL, filename = NULL, path = ".", signatures = NULL, saveHTML = FALSE, filename.out = NULL, X11 = TRUE, verbose = FALSE)
}

\arguments{
  \item{data}{a DBFMCLresult class object or a matrix.}
  \item{filename}{a character string representing the file name (if data are loaded from a flat file).}
  \item{path}{a character string representing the data directory.}
  \item{signatures}{a vector that indicates which signatures to plot.}
  \item{saveHTML}{if set to TRUE a HTML file is created.}
  \item{filename.out}{a character string representing the HTML file name (if \code{saveHTML=TRUE}).}
  \item{X11}{if set to TRUE a new graphic windows is generated.}
  \item{verbose}{if set to TRUE the function runs verbosely.}
}
\details{
Mean expression profile is highlighted in green.
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\seealso{\code{\link{matplotProfiles}}}
\examples{
\dontrun{
# Please check that the directory is writable
# before using the following code 
library(ALL)
data(ALL)
ALLnorm <- doNormalScore(exprs(ALL))
res <- DBFMCL(data=ALLnorm, name="ALLout")
plotGeneExpProfiles(res)
plotGeneExpProfiles(res, signatures=1:2)
plotGeneExpProfiles(res, signatures=1:2, saveHTML=TRUE)
}
}

\keyword{hplot}
