\name{getData4DBFMCL}
\alias{getData4DBFMCL}

\title{Fetch an expression matrix from a file or a Seurat object}
\description{This function retrieves the expression matrix from a file or a \code{Seurat} object. This is an internal function called by createSignatures4TB that should not be used directly.}
\usage{
getData4DBFMCL(data = NULL, filename = NULL, path = ".")
}

\arguments{
  \item{data}{a matrix}
  \item{filename}{a character string representing the file name.}
  \item{path}{A character string representing the data directory.}
}

\section{Warnings}{
Convert data.frame, Seurat object or tab-delimited file to matrix class object. The input data must contain an expression matrix with gene as rows and samples as columns. Note that space characters inside gene names are not allowed (as they are not supported by the mcl command-line program).
}

\value{
This function return a list which contains a matrix and a name. If filename is not equal to NULL, this name corresponds to the prefix of the given filename, in the other case name = NULL.
}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}
\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}


\seealso{\code{\link{createSignatures4TB}}}


\keyword{manip}