\name{heatmapFromCDT}
\alias{heatmapFromCDT}

\title{Reads a "cdt" file containing a set of transcriptional signatures and produces heatmaps. }
\description{This function can be called to display the results of the createSignatures4TB function.}

\usage{
heatmapFromCDT(cdt.filename, signature = NULL, fac = NULL)
}

\arguments{
  \item{cdt.filename}{a character string representing the "cdt" file name.}
  \item{signature}{a vector. By default (\code{signature = NULL}), all signatures are displayed.}
  \item{fac}{a factor. By default \code{fac = NULL}, no phenodata are displayed. If fac corresponds to a factor, the several levels of fac are displayed by a distinct color.}
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\seealso{\code{\link{createSignatures4TB}}, \code{\link{plotGeneExpProfiles}}}

\examples{
\dontrun{
library(ALL)
data(ALL)
res <- createSignatures4TB(data=ALL, name="ALLdataset", median.center=TRUE, distance.method="pearson")
all <- heatmapFromCDT("ALLdataset.dataMods.cdt")
sig6 <- heatmapFromCDT("ALLdataset.dataMods.cdt", signature=6)
rownames(sig6)
}}
\keyword{hplot}
