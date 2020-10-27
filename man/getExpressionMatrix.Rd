\name{getExpressionMatrix}
\alias{getExpressionMatrix}

\title{A function to fetch an expression matrix from the TranscriptomeBrowser database.}
\description{
This function takes a signature ID as input and retrieves the corresponding expression matrix from the TranscriptomeBrowser database (\url{http://tagc.univ-mrs.fr/tbrowser}).}

\usage{
getExpressionMatrix(signatureID = NULL, verbose = TRUE, save = FALSE)
}

\arguments{
  \item{signatureID}{a TBrowserDB signature ID (e.g., "0AC1A39A5"), which can be obtained using the \code{\link{getSignatures}} function.}
  \item{verbose}{if set to TRUE the function runs verbosely.}
  \item{save}{if set to TRUE data are stored onto disk.}
}

\value{
A \code{data.frame} containing the expression matrix of the requested signature. The first two columns store probe and gene informations and additional columns contain expression values. for corresponding samples (as columns) 
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\examples{
\dontrun{
## what are the TBrowserDB signatures that contain both XBP1, GATA3 and ESR1 ?
library(RColorBrewer)
res <- getSignatures(field="gene", value="XBP1 & GATA3 & ESR1")
res

# Retrieve expression matrix for the second signature.
em <- getExpressionMatrix(signatureID = "3DE64836D")

#Getting gene names and sample informations
geneNames <- paste(em[,1],em[,2],sep="||") 
em <- as.matrix(em[,-c(1,2)])
samplesInfo <- getTBInfo(field="samples", value = "3DE64836D") 
ind <- match(colnames(em), samplesInfo[,1])
colnames(em) <- samplesInfo[ind,2]

# Heatmap
col <- colorRampPalette(brewer.pal(10, "RdBu"))(256) 
split <- strsplit(colnames(em)," (", fixed=TRUE) 
pheno <- unlist(lapply(split,"[",1))
pheno <- as.factor(pheno)
levels(pheno) <- 1:5
cc <- rainbow(5, start=0, end=.3)
cc <- cc[pheno]
heatmap(em, Rowv=FALSE, Colv=FALSE, col=col, ColSideColors=cc, labRow=geneNames,cexRow=0.3)
}
}
\keyword{manip}
