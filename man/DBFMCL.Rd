\name{DBFMCL}
\alias{DBFMCL}

\title{The "Density Based Filtering and Markov CLustering" algorithm (DBF-MCL).}
\description{DBF-MCL is a tree-steps adaptative algorithm (\url{http://tagc.univ-mrs.fr/tbrowser/}) that \emph{(i)}find elements located in dense areas (DBF), \emph{(ii)}uses selected items to construct a graph, \emph{(iii)}performs graph partitioning using the Markov CLustering Algorithm (MCL). 

This function requires installation of the mcl program (\url{http://www.micans.org/mcl}). See "Warnings" section for more informations.
}


\usage{
DBFMCL(data = NULL, filename = NULL, path = ".", name = NULL, distance.method = c("pearson", "spearman", "euclidean", "spm", "spgm"), 
 clustering = TRUE, silent = FALSE, verbose = TRUE, k = 150, random = 3, memory.used = 1024, fdr = 10, inflation = 2.0, set.seed = 123, returnRank = FALSE)
}

\arguments{
  \item{data}{a \code{matrix}, \code{data.frame} or \code{ExpressionSet} object.}
  \item{filename}{a character string representing the file name.}
  \item{name}{a prefix for the names of the intermediary files created by DBF and MCL.}
  \item{path}{a character string representing the data directory where intermediary files are to be stored. Default to current working directory.}
  \item{distance.method}{a method to compute the distance to the k-th nearest neighbor. One of "pearson" (Pearson's correlation coefficient-based distance), "spearman" (Spearman's rho-based distance), "euclidean". The "spm" distance corresponds to the arithmetic mean :("pearson"+"spearman")/2 whereas "spgm" is the geometric mean : sqrt("pearson"*"spearman).}
  \item{clustering}{indicates whether partitioning step (MCL) should be applied to the data. If \code{clustering = FALSE}, the function returns a \code{DBFMCLresult} object that contains informative elements (as detected by the DBF step) coerced into a single cluster.   }
  \item{silent}{if set to TRUE, the progression of distance matrix calculation is not displayed.}
  \item{verbose}{if set to TRUE the function runs verbosely.}
  \item{k}{the neighborhood size.}
  \item{random}{the number of simulated distributions S to compute. By default \code{random = 3}.}
  \item{memory.used}{size of the memory used to store part of the distance matrix. The subsequent sub-matrix is used to computed simulated distances to the k-th nearest neighbor (see detail section).}
  \item{fdr}{an integer value corresponding to the false discovery rate (range: 0 to 100).}
  \item{inflation}{the main control of MCL. Inflation affects cluster granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}. \code{inflation = 5.0} will tend to result in fine-grained clusterings, and whereas \code{inflation = 1.2} will tend to result in very coarse grained clusterings. By default, \code{inflation = 2.0}. Default setting gives very good results for microarray data when k is set between 70 and 250.}
  \item{set.seed}{specify seeds for random number generator.}
  \item{returnRank}{allows to obtain in the \code{DBFMCLresult} object, a rank's matrix. The output files are processed using the \code{normalization} argument.}
}
\details{
When analyzing a noisy dataset, one is interested in isolating dense regions as they are populated with genes/elements that display weak distances to their nearest neighbors (i.e. strong profile similarities). To isolate these regions DBF-MCL computes, for each gene/element, the distance with its kth nearest neighbor (DKNN).In order to define a critical DKNN value that will depend on the dataset and below which a gene/element will be considered as falling in a dense area, DBF-MCL computes simulated DKNN values by using an empirical randomization procedure. Given a dataset containing n genes and p samples, a simulated DKNN value is obtained by sampling n distance values from the gene-gene distance matrix D and by extracting the kth-smallest value. This procedure is repeated n times to obtain a set of simulated DKNN values S. Computed distributions of simulated DKNN are used to compute a FDR value for each observed DKNN value. The critical value of DKNN is the one for which a user-defined FDR value (typically 10\%) is observed. Genes with DKNN value below this threshold are selected and used to construct a graph. In this graph, edges are constructed between two genes (nodes) if one of them belongs to the k-nearest neighbors of the other. Edges are weighted based on the respective coefficient of correlation (\emph{i.e.}, similarity) and the graph obtained is partitioned using the Markov CLustering Algorithm (MCL). 
}

\value{
a DBFMCLresults class object.
}

\section{Warnings}{
With the current implementation, this function only works only on UNIX-like plateforms.

MCL should be installed. One can used the following command lines in a terminal:

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
}

\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}
\seealso{\code{\link{createSignatures4TB}}}
\examples{
\dontrun{
## with an artificial dataset

m <- matrix(rnorm(80000), nc=20)
m[1:100,1:10] <- m[1:100,1:10] + 4
m[101:200,11:20] <- m[101:200,11:20] + 3
m[201:300,5:15] <- m[201:300,5:15] + -2

res <- DBFMCL(data = m, distance.method = "pearson", clustering = TRUE, k = 25)
plotGeneExpProfiles(res)
plotGeneExpProfiles(res,signatures=1)

## with a real dataset
library(ALL)
data(ALL)
sub <- exprs(ALL)[1:3000,]

#First, we will normalize the data set using the doNormalScore function. 
subNorm <- doNormalScore(sub)
res <- DBFMCL(subNorm, distance.method="pearson", memory=512)

#The results are stored in an instance of class DBFMCLresult.
class(res)
res

# The expression matrix is stored in the data slot. 
# This matrix contains only genes detected as informative (that is falling into a cluster).
head(res@data[,1:2])

# The partitioning results are stored in the cluster slot.
slotNames(res)

# Here, 3 TS were found.
res@size

# The following instruction can be used to get the expression matrix corresponding to the first TS. 
res@data[res@cluster ==1,]

# The high level function plotGeneExpProfilescan be used to visualize, 
# for instance, gene expression profiles corresponding to the first signature.
plotGeneExpProfiles(res, sign=1)

#To stored the partitioning results onto disk (as a tab-delimited file), 
# use the writeDBFMCLresult function as show below.
writeDBFMCLresult(res, filename.out="ALL.sign.txt")

}
}

\keyword{manip}

