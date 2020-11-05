\name{getTBInfo}
\alias{getTBInfo}
\title{Get information about an experiment or a platform.}
\description{This function fetch informations from the TBrowserDB for a microarray experiment, a microarray platform or a transcriptional signature.}
\usage{
getTBInfo(field = c("platform", "experiment", "signature", "samples"), value = NULL, verbose = TRUE, save = FALSE)
}

\arguments{
  \item{field}{The type of information to retrieve. Should be one of "platform", "experiment" or "signature".}
  \item{value}{A platform ID (e.g., GPL96), experiment ID (e.g., GSE2004) or signatureID (e.g., "0AC1A39A5") depending on the "field" argument.}
  \item{verbose}{If set to TRUE the fonction displays informations on the screen.}
  \item{save}{if set to TRUE data are stored onto disk.} 
}

\value{
The output will differ depending one the \code{field} argument.
  \item{}{ if \code{field = "platform"}, the function will display plateform informations (e.g.: name, organism, manufacturer, nb. probes, nb. genes, Title and Description).}
  \item{}{ if \code{field = "experiment"}, the function will display informations about the experiment (e.g.: name, organism, PMID, nb. samples, title and summary).}
  \item{}{ if \code{field = "signature"}, the function will display informations about the experiment (e.g.: signatureID, platform, experiment, organism, nb.probes, nb.genes, nb.samples).}
  \item{}{ if \code{field = "samples"}, the function will display informations about all the samples from an signatureID (e.g.: sampleID, descriptions).}
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\seealso{Other function which allow to query the TBrowser database: \code{\link{getSignatures}}, \code{\link{getExpressionMatrix}}}
\examples{
\dontrun{
# retrieving information related to GSE2004
gse2004Info <- getTBInfo(field="experiment", value="GSE2004")

# retrieving information related to GPL96
gpl96Info <- getTBInfo(field="platform", value="GPL96")

# retrieving information related to a signature
signInfo <- getTBInfo(field="signature", value="0AC1A39A5")

# retrieving samples information related to a signature
samplesInfo <- getTBInfo(field="samples", value="0AC1A39A5")
}
}

\keyword{manip}
