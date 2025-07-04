% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{discrete_palette}
\alias{discrete_palette}
\title{Generate a discrete color palette}
\usage{
discrete_palette(n = 10, palette = c("Ju1", "De1", "ggplot"))
}
\arguments{
\item{n}{An integer specifying the number of colors to generate.}

\item{palette}{A character vector specifying the palette to use.}
}
\value{
A character vector of color codes.
}
\description{
This function generates a vector of colors for a discrete variable, 
given a specified palette name.
}
\examples{
discrete_palette()
discrete_palette(n=20, palette = "ggplot")
}
