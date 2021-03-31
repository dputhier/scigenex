\name{writeDBFMCLresult}
\alias{writeDBFMCLresult}

\title{Writes the results of DBFMCL function onto disk.}
\description{Writes the results of DBFMCL function onto disk. The output file contains an expression matrix in which transcriptional signatures are separed by blank lines.}

\usage{
writeDBFMCLresult(object, filename.out = NULL, path = ".", verbose = TRUE)
}

\arguments{
  \item{object}{a \code{DBFMCLresults} class object. }
  \item{filename.out}{a character string representing the file name where the data will be stored.}
  \item{path}{a character string representing the data directory where the data will be stored.}
  \item{verbose}{if set to TRUE the function runs verbosely.}
}
\value{
A tab-delimited file.
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\seealso{\code{\link{DBFMCLresult-class}}, \code{\link{DBFMCL}} }
\examples{
\dontrun{
m <- matrix(rnorm(80000), nc=20)
m[1:100,1:10] <- m[1:100,1:10] + 4
m[101:200,11:20] <- m[101:200,11:20] + 3
m[201:300,5:15] <- m[201:300,5:15] + -2
writeDBFMCLresult(res, "ALL.sign.txt")
}



}

\keyword{manip}