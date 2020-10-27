\name{RTools4TB_1.1.4-package}
\alias{RTools4TB_1.1.4-package}
\alias{RTools4TB}
\docType{package}
\title{The RTools4TB package: data mining of public microarray data through connections to the TranscriptomeBrowser database.}

\description{
TranscriptomeBrowser (TBrowser, \url{http://tagc.univ-mrs.fr/tbrowser}) hosts a large collection of transcriptional signatures (TS)  automatically extracted from the Gene Expression Omnibus (GEO) database. Each GEO experiment (GSE) was processed so that a subset of the original expression matrix containing the most relevant/informative genes was kept and organized into a set of homogeneous signatures. Each signature was tested for functional  enrichment using annotations terms obtained from numerous ontologies or curated databases (Gene Ontology, KEGG, BioCarta, Swiss-Prot, BBID, SMART, NIH Genetic Association DB, COG/KOG...) using the DAVID knowledgebase. The RTools4TB package can be used to perform complexe queries to the database. RTools4TB can be helpful (i) to define the biological contexts (i.e, experiments) in which a set of genes are co-expressed and (ii) to define their most frequent neighbors.

In addition, RTools4TB comes with a new algoritm, "Density Based Filtering And Markov Clustering" (DBF-MCL), whose goal is to partition large and noisy datasets. DBF-MCL is a tree-step adaptative algorithm that (i) find elements located in dense areas (ie. clusters) (ii) uses selected items to construct a graph and (iii) performs graph partitioning using MCL. This algorithm is implemented in the RTools4TB package although it requires a UNIX-like systems.
}

\details{
\tabular{ll}{
Package: \tab RTools4TB\cr
Type: \tab Package\cr
Version: \tab 1.1.43\cr
Date: \tab 2008-10-08\cr
License: \tab LGPL\cr
Notes: \tab Please note that a web connection is required to fetch expression data from the TrancriptomeBrowser database.\cr
}
}
\author{
Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.
Maintainer: Aurelie Bergon \email{bergon@tagc.univ-mrs.fr}
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}

\keyword{package}
