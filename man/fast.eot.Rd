% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fast.eot.R
\name{fast.eot}
\alias{fast.eot}
\title{A parallelized version to calculate  EOTs (Empirical orthogonal teleconnections)}
\usage{
fast.eot(
  x,
  n = 1,
  standardised,
  n_cores = parallel::detectCores() - 1,
  verbose = T
)
}
\arguments{
\item{x}{The spatio-temporal data, organised in a 'RasterBrick'-format.}

\item{n}{The number of EOTs to calculate.}

\item{standardised}{Boolean. Indicates whether the output should be scaled or not.}

\item{n_cores}{Integer. The number of CPU-cores for the function to use.}

\item{verbose}{Boolean. Should intermediate messages be printed to the console (TRUE) or not (FALSE). Default to TRUE.}
}
\value{
A List of the EOTs.
}
\description{
This function builds on the work done for the 'remote' package, but parallelizes the code, in order to reduce computation times.
Given a spatio-temporal field (as a raster/brick file) it calculates a predefined number of EOTs.
}
\references{
\bold{Empirical Orthogonal Teleconnections}\cr
H. M. van den Dool, S. Saha, A. Johansson (2000)\cr
Journal of Climate, Volume 13, Issue 8, pp. 1421-1435\cr
\url{http://journals.ametsoc.org/doi/abs/10.1175/1520-0442%282000%29013%3C1421%3AEOT%3E2.0.CO%3B2
}

\bold{Empirical methods in short-term climate prediction}\cr
H. M. van den Dool (2007)\cr
Oxford University Press, Oxford, New York\cr
\url{https://global.oup.com/academic/product/empirical-methods-in-short-term-climate-prediction-9780199202782?cc=de&lang=en&}
}
