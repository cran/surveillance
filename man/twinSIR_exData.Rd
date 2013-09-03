\encoding{latin1}
\name{twinSIR_exData}
\alias{foodata}
\alias{fooepidata}
\alias{foofit}
%\alias{bw}
\docType{data}

\title{
Artificial data and data from the German Federal State Baden-Wuerttemberg
}

\description{
  This is data, which is used in \code{twinSIR} examples or the package
  vignette.
  
  \code{fooepidata} is an object of class \code{"\link{epidata}"} simulated from
  artificial \code{foodata}. \code{foofit} is an object of class
  \code{"\link{twinSIR}"} and is the result of a fit of \code{fooepidata}.
  
  % The data set \code{bw} contains some information on the counties of
  % Baden-\enc{Württemberg}{Wuerttemberg} for the years 2002--2006:
  % latitude and longitude of the county's centroid, population, area,
  % annual number of infections with tick-borne encephalitis (where the county of
  % infection is known), as well as some yearly aggregated weather information
  % (temperature, precipitation, sunshine).
}

\examples{
data(foodata)
data(fooepidata)
data(foofit)
}

\keyword{datasets}
\keyword{internal}
