\name{colorScale}
\alias{colorScale}

\title{Calibration bar for color images}
\description{This function is the same as the \code{\link[marray]{maPalette}} function in the marray library. It is used to produce a color scale.
}
\usage{
colorScale(low = "white", high = c("green", "red"), mid = NULL, k = 50)
}

\arguments{

\item{low}{Color for the lower end of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of colors), a hexadecimal string of the form \code{#rrggbb}, or an integer i meaning \code{palette()[i]}}
\item{high}{Color for the upper end of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of colors), a hexadecimal string of the form \code{#rrggbb}, or an integer i meaning \code{palette()[i]}}
\item{mid}{Color for the middle portion of the color palette, specified using any of the three kinds of R colors, i.e., either a color name (an element of colors), a hexadecimal string of the form \code{#rrggbb}, or an integer i meaning \code{palette()[i]}}
\item{k}{Number of colors in the palette}
}

\author{see the \code{\link[marray]{maPalette}} }

\examples{
Rcol <- colorScale(low="white", high="red", k=10)
Gcol <- colorScale(low="white", high="green", k=50)
RGcol <- colorScale(low="green", high="red", k=100)
}


\keyword{aplot}
\keyword{hplot}
