\name{MCL}
\alias{MCL}

\title{Invokes the Markov CLustering algorithm (MCL).}
\description{This function invokes the mcl system command. MCL is a clustering algorithm for graphs that was developped by Stijn van Dongen (see references for further informations).}
\usage{
MCL(name, inflation = 2.0, silent = FALSE)
}

\arguments{
  \item{name}{a character string corresponding to the file name.}
  \item{inflation}{the main control of MCL. Inflation affects cluster granularity. It is usually chosen somewhere in the range \code{[1.2-5.0]}. \code{inflation = 5.0} will tend to result in fine-grained clusterings, and whereas \code{inflation = 1.2} will tend to result in very coarse grained clusterings. By default, \code{inflation = 2.0}. Default setting gives very good results for microarray data when k is set around 100.}
  \item{silent}{if set to TRUE, the progression of the MCL partitionning is not displayed.}
}

\value{
Returns a file with the ".mcl\_out.txt" extension.
}

\section{warning}{Works only on UNIX-like plateforms. MCL should be installed. The following command lines can be used for installation.

\code{# Download the latest version of mcl (RTools4TB has been tested successfully with the 06-058 version).}

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
Stijn  van Dongen. A cluster algorithm for graphs.  Technical Report INS-R0010, National Research Institute for Mathematics and Computer Science in the Netherlands, Amsterdam, May 2000. \url{http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z}

Lopez F.,Textoris J., Bergon A., Didier G., Remy E., Granjeaud S., Imbert J. , Nguyen C. and Puthier D. TranscriptomeBrowser: a powerful and flexible toolbox to explore productively the transcriptional landscape of the Gene Expression Omnibus database. PLoSONE, 2008;3(12):e4001.
}
\author{Bergon A., Lopez F., Textoris J., Granjeaud S. and Puthier D.}

\keyword{manip}

