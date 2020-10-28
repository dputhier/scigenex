\name{createSignatures4TB}
\alias{createSignatures4TB}

\title{Creates a set of transcriptional signatures from a microarray dataset.}
\description{This function is a wrapper to create sets of transcriptional signatures (as in the TranscriptomeBrowser Project, TBrowser,\url{http://tagc.univ-mrs.fr/tbrowser}). This function creates a "cdt" file containing a set of expression matrices (transcriptional signatures) separated by blank lines. Please note that it requires both MCL and Cluster 3.0 (see 'warnings section'). It accepts both a matrix or file name as input.
}


\usage{
createSignatures4TB(data = NULL, filename = NULL, path = ".", name = NULL, normalization = c("rank", "gaussian", "quantiles", "none"),
 distance.method = c("spearman", "pearson", "euclidean", "spm", "spgm"), silent = FALSE, verbose = TRUE, k = 150, random = 3, memory.used = 1024, fdr = 10,
 inflation = 2.0, median.center = FALSE, set.seed = 123, returnRank = FALSE)
}

\arguments{
  \item{data}{a \code{matrix}, \code{data.frame} or \code{ExpressionSet} object.}
  \item{filename}{a character string representing the file name to load.}
  \item{path}{a character string representing the data directory.}
  \item{name}{a prefix for the name of the created files.}
  \item{normalization}{indicates whether data should be normalized prior to analysis (see details).}
  \item{distance.method}{a method to compute the distance to the k-th nearest neighbor. One of "pearson" (Pearson's correlation coefficient-based distance), "spearman" (Spearman's rho-based distance), "euclidean", "spm" or "spgm". Note that the "spm" distance corresponds to the arithmetic mean of pearson- and spearman-based distance : ("pearson"+"spearman")/2 whereas "spgm" computes their geometric mean : sqrt("pearson"*"spearman").}
  \item{silent}{if set to TRUE, the progression of distance matrix calculation is not displayed.}
  \item{verbose}{if set to TRUE the function runs verbosely.}
  \item{k}{the neighborhood size.}
  \item{random}{the number of simulated distributions S to compute. By default \code{random = 3}.}
  \item{memory.used}{size of the memory used to store part of the distance matrix. The subsequent sub-matrix is used to computed simulated distances to the k-th nearest neighbor (see detail section of \code{\link{DBFMCL}} function).}
  \item{fdr}{an integer value corresponding to the false discovery rate (range: 0 to 100).}
  \item{inflation}{the main control of MCL. Inflation affects cluster granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}. \code{inflation = 5.0} will tend to result in fine-grained clusterings whereas \code{inflation = 1.2} will tend to result in very coarse grained clusterings. By default, \code{inflation = 2.0}. Default setting gives very good results for microarray data.}
  \item{median.center}{if set to TRUE, median-centering is applied to the rows of the matrix.}
  \item{set.seed}{specify seeds for random number generator.}
  \item{returnRank}{This argument modifies the output. Given a set of elements conserved after the filtering step of the DBFMCL algorithm, if \code{returnRank = TRUE} their expression values are replaced by their corresponding ranks in the input matrix.}
}

\details{
The Markov Cluster Algorithm was written by S. Van Dongen (see reference section). Cluster was originally written by Michael Eisen (\url{http://rana.lbl.gov/EisenSoftware.htm}). The command line version of Cluster version 3.0 was created by Michiel de Hoon, together with Seiya Imoto and Satoru Miyano. 
}

\section{Warnings}{
With the current implementation, this function works only on UNIX-like plateforms.

Cluster 3.0 should be installed in its command-line only version:

Please see \url{http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/software.htm} for further informations.

\code{wget http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/cluster-1.36.tar.gz}

\code{tar xvfz cluster-1.36.tar.gz}

\code{cd cluster-1.36/}

\code{./configure --without-x}

\code{make}

\code{sudo make install}

\code{# You should get cluster in your path}

\code{cluster -v}


MCL should be installed:

\code{# Download the latest version of mcl (the script has been tested successfully with the 06-058 version).}

\code{wget http://micans.org/mcl/src/mcl-latest.tar.gz}

\code{# Uncompress and install mcl}

\code{tar xvfz mcl-latest.tar.gz}

\code{cd mcl-xx-xxx}

\code{./configure}

\code{make}

\code{sudo make install}

\code{# You should get mcl in your path}

\code{mcl -h}


}

\references{
Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and 
flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.

Van Dongen S. (2000) A cluster algorithm for graphs. National Research Institute for Mathematics and Computer Science in the 1386-3681.

Open source clustering software. De Hoon MJ, Imoto S, Nolan J, Miyano S. Bioinformatics. 2004 Jun 12;20(9):1453-4.

Cluster analysis and display of genome-wide expression patterns. Eisen MB, Spellman PT, Brown PO, Botstein D. Proc Natl Acad Sci U S A. 1998 Dec 8;95(25):14863-8.
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}
\seealso{\code{\link{DBFMCL}},\code{\link{heatmapFromCDT}},\code{\link{plotGeneExpProfiles}},\code{\link{getSignatures}},\code{\link{getExpressionMatrix}}}

\examples{
\dontrun{
## with an artificial dataset

m <- matrix(rnorm(80000), nc=20)
m[1:100,1:10] <- m[1:100,1:10] + 4
m[101:200,11:20] <- m[101:200,11:20] + 3
m[201:300,5:15] <- m[201:300,5:15] + -2

res <- createSignatures4TB(data = m, name="artificial", distance.method = "pearson", median.center=TRUE, k = 25)

plotGeneExpProfiles(res)
allsign <- heatmapFromCDT("artificial.dataMods.cdt")

plotGeneExpProfiles(res, signature=1)
heatmapFromCDT("artificial.dataMods.cdt", signature=1)

## with a real dataset
library(ALL)
data(ALL)
exp <- createSignatures4TB(data = ALL , name="ALLdataset", distance.method = "pearson", median.center=TRUE, k = 100)
plotGeneExpProfiles(exp, signatures=1)
plotGeneExpProfiles(res)
allsign <- heatmapFromCDT("ALLdataset.dataMods.cdt")
si1 <- heatmapFromCDT("ALLdataset.dataMods.cdt", signature=1)
}
}

\keyword{manip}
