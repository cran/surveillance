\name{measels.weser}
\alias{measels.weser}
\docType{data}
\title{Measels epidemics in Lower Saxony in 2001-2002}
\description{
  Weekly counts of new measels cases for each "Kreis" of the 
  administrative district "Weser-Ems" in Lower Saxony, Germany, in 2001 and 2002.
  All in all there are 15 "Kreise", two "Kreise" have been omitted %due to ?
  
  }
\usage{data(measels.weser)}
\format{
  An multivariate object of class disProg with 104 observations for each one of the 15 Kreise.
  \describe{
    \item{week}{Number of week.}
    \item{observed}{Matrix with number of counts in the corresponding week and Kreis.}
    \item{state}{Boolean whether there was an outbreak -- dummy not implemented.}
    \item{neighbourhood}{Neighbourhood matrix.}
    \item{populationFrac}{Population fractions.} %
  }
}
\source{
  ??
}
\examples{
data(measels.weser)
plot(measels.weser, as.one=FALSE) 
}
\keyword{datasets}
