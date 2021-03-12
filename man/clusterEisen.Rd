\name{clusterEisen}
\alias{clusterEisen}

\title{Performs system call to the Cluster 3.0 program.}
\description{This function does a system call to the Cluster 3.0 program. It requires Cluster 3.0 to be installed on the machine and to be available in the path (as "cluster"). This internal function is used by \code{\link{createSignatures4TB}}. it performs only clustering of the columns (using a Pearson's correlation coefficient-based distance). 
}
\usage{
clusterEisen(filename, median.center = FALSE, silent = FALSE)
}

\arguments{
  \item{filename}{a character string representing the file name.}
  \item{median.center}{a logical indicating whether rows should be median-centered.}
  \item{silent}{if set to TRUE, the progression of clustering is not displayed.}
}
\details{
Cluster was originally written by Michael Eisen. The command line version of Cluster version 3.0 (for Windows, Mac OS X, Unix, and Linux) was created by Michiel de Hoon, together with Seiya Imoto and Satoru Miyano.
}

\section{Warnings}{
Only tested on UNIX-alikes plateforms.

Cluster 3.0 should be installed in its command-line only version:

\code{wget http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/cluster-1.36.tar.gz}

\code{tar xvfz cluster-1.36.tar.gz}

\code{cd cluster-1.36/}

\code{./configure --without-x}

\code{make}

\code{sudo make install}

\code{# You should get cluster in your path}

\code{cluster -v}


Please see \url{http://bonsai.ims.u-tokyo.ac.jp/~mdehoon/software/cluster/software.htm} for further informations.

}

\references{
Open source clustering software. De Hoon MJ, Imoto S, Nolan J, Miyano S. Bioinformatics. 2004 Jun 12;20(9):1453-4.

Cluster analysis and display of genome-wide expression patterns. Eisen MB, Spellman PT, Brown PO, Botstein D. Proc Natl Acad Sci U S A. 1998 Dec 8;95(25):14863-8.
}


\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}


\keyword{manip}
