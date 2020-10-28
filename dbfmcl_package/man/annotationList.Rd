\name{annotationList}
\alias{annotationList}
\docType{data}
\title{The keywords used to get information about functional enrichment of transcriptional signatures.}
\description{A dataframe containing the keywords used to get information about functional enrichment of transcriptional signatures.}
\usage{data(annotationList)}
\examples{
data(annotationList)
names(annotationList)
attach(annotationList)
table(TableName)
annotationList[1:4,]
head(annotationList[TableName=="KEGG_PATHWAY",])
}
\keyword{datasets}
