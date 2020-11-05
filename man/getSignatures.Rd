\name{getSignatures}
\alias{getSignatures}
\title{A function to retrieve transcriptional signature IDs from the TranscriptomeBrowser database (TBrowserDB).}

\description{
This is one of the main function of the \pkg{RTools4TB} package. It allows direct access to TBrowserDB (\url{http://tagc.univ-mrs.fr/tbrowser}). The \code{getSignatures} function can be used to retrieve transcriptional signatures (i)derived from a given experiment or microarray platform, (ii)containing a user-defined list of genes or probes (using or not a boolean query) or (iii)enriched in genes sharing a common annotation term (user must provide a q-value).

See "Details" section for more information about the syntax.}

\usage{
getSignatures(field=c("gene", "probe", "platform", "experiment", "annotation"), value = NULL, qValue = NULL, nbMin = NULL, verbose = TRUE, save = FALSE)
}

\arguments{

    \item{field}{The request type. Should be one of: "gene", "probe", "platform", "experiment", "annotation".}
    \item{value}{Depends on the \code{field} argument:
    if \code{field} is set to "gene" \code{value} must contain HUGO IDs (\emph{e.g.}, \code{"CD4 CD3E CD3D"}). Logical operators are supported (\emph{e.g.}, \code{"CD4 & CD3E & CD3D"}, see "details" section).
    if \code{field} is set to "probe" \code{value} must contain a list of probe IDs (\emph{e.g.}, Affymetrix probe IDs). Logical operators are supported.
    if \code{field} is set to "platform" \code{value} must contain one platform ID (\emph{e.g.}, \code{"GPL96"}).
    if \code{field} is set to "experiment" \code{value} must contain one experiment ID (\emph{e.g.}, \code{"GSE2004"}).
    if \code{field} is set to "annotation" \code{value} must contain a list of annotation terms separated by logical operators(\emph{e.g.}, \code{"breast cancer"} or \code{"18q11.2|18q12.1|18q21.1|18q22-q23"}). 
   }
  \item{qValue}{an integer (10E-"qValue"). Default to 0. This q-value is used to select signatures associated with a given annotation term (see examples section). Used only when \code{field = "annotation"}.}
   \item{nbMin}{an integer. Used only when \code{value} corresponds to a gene list without logical operators (see details). Only signatures containing at least nbMin genes out of the list will be retrieved (see details section).}
  \item{verbose}{if set to TRUE the function runs verbosely.}
  \item{save}{if set to TRUE data are stored onto disk.}  
}
\details{
The "value" argument to getSignatures may contain logical operators (see help section on TBrowser web site for more informations, \url{http://tagc.univ-mrs.fr/tbrowser})

\code{&}    : AND
\code{|}    : OR
\code{!}    : NOT , (used in conjonction with &)

However, when \code{field = "gene"} or \code{field = "probe"}, user can perform a request using a list of item separated by blanks (without logical operators). These blanks are interpreted as the OR logical operators. In this case, all signatures containing at least one gene of the list will be returned. To select more informative signatures we suggest to use the \code{nbMin} argument that will select signatures containing at least \code{nbMin} genes out of the list.

Moreover, user may include logical operators in the request. Indeed, this is a convenient way to create relevant queries. Suppose your field of interest is related to T-cell activation. You could be interested in retrieving all TS that contain the CD4 gene as they should contain additional T cell markers. Comparing these TS should help you to define a set of frequent CD4 neighbors (very likely related to TCR signaling cascade). Thereby, your request should be:

\code{res <- getSignatures(field="gene", value="CD4")}

This gene is found in 371 TS (with the current database release), and obtaining associated gene lists would be time consuming and would not emphasize on what you are really expecting. Indeed, the CD4 marker is also expressed by macrophages. Another solution would be to search for TS containing two T-cell markers (CD4 and CD3E for instance) and to exclude (using the NOT operator) those containing the CD14 marker (a macrophages marker). The syntax should be the following:

\code{res <- getSignatures(field="gene", value="CD4 & CD3E & !CD14")}

In the same way you could try to exclude TS containing B-cells by discarding those containing the CD19 of IGHM marker. The resulting query would be the following:

\code{res <- getSignatures(field="gene", value="CD4 & CD3E & !(CD19 | IGHM)")}

}
\value{
This function will return a vector containing the names of the transcriptional signatures that satisfy the constraints. Additional informations about these signatures (GEO platform ID, GEO experiment ID, Organism, number of probes, number of genes, number of biological samples) can be obtained using the \code{\link{getTBInfo}} function (\code{field = "signatureID"}). 
}
\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}
\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\seealso{Other functions which allow to query the TBrowser database: \code{\link{getTBInfo}}, \code{\link{getExpressionMatrix}}}
\examples{
\dontrun{
# retrieving transcriptional signatures containing PCNA, CDC2 and CDC6.
res <- getSignatures(field="gene", value="PCNA & CDC2 & CDC6")

# retrieving transcriptional signatures contain at least two genes out of the following list: "PCNA, CDC2 and CDC6".
res <- getSignatures(field="gene", value="PCNA CDC2 CDC6", nbMin=2)

# retrieving transcriptional signatures related to GSE2004
gse2004TS <- getSignatures(field="experiment", value="GSE2004")

# retrieving transcriptional signatures related to the platform GPL96
gpl96TS <- getSignatures(field="platform", value="GPL96")

# retrieving transcriptional signatures enriched in gene related to the keyword ""HSA04110:CELL CYCLE" (KEGG_PATHWAY)
data(annotationList)
attach(annotationList)
table(TableName)
annotationList[Keyword=="HSA04110:CELL CYCLE",]
ccTS20 <- getSignatures(field="annotation", value="HSA04110:CELL CYCLE", qValue=20)


# retrieving transcriptional signatures enriched in gene located in 8q region.
query <- paste(grep("^8q", Keyword, val = T), collapse = "|")
query
cc <- getSignatures(field = "annotation", value = query, qValue = 10)

}
}
\keyword{manip}
