\name{momo}
\alias{momo}
\docType{data}
\encoding{latin1}

\title{Danish 1994-2008 all cause mortality data for six age groups}
\description{
  Weekly number of all cause mortality from 1994-2008 in each
  of the six age groups <1, 1-4, 5-14, 15-44, 45-64, 65-74,
  75-84 and 85 years. 
}
\usage{data(momo)}
\details{
  The object of class \code{"\linkS4class{sts}"} contains the number of all cause
  mortality from 1994-2008 in Denmark for each
  of the six age groups <1, 1-4, 5-14, 15-44, 45-64, 65-74,
  75-84 and 85 years. A special feature of such EuroMOMO data is that
  weeks are handled as defined by the ISO 8601 standard, which can be
  handled by the \code{"sts"} class.

  The \code{population} slot of the \code{momo} object contains the
  population size in each of the six age groups.
  These are yearly data obtained from the StatBank Denmark.

  The aim of the EuroMOMO project is to develop and strengthen real-time
  monitoring of mortality across Europe; this will enhance the
  management of serious public health risks such as pandemic influenza,
  heat waves and cold snaps. For further details see the homepage of the
  EuroMOMO project. 
}
\source{
  Department of Epidemiology, Statens Serum Institute, Copenhagen, Denmark
  StatBank Denmark, Statistics Denmark, \url{http://www.statistikbanken.dk/}
}
\examples{
data("momo")
plot(momo)
}
\references{
  H\enc{ö}{oe}hle, M. and A. Mazick, A. (2009) Aberration detection in R
  illustrated by Danish mortality monitoring,  Book chapter to appear in
  T. Kass-Hout and X. Zhang (Eds.) Biosurveillance: A Health Protection
  Priority, CRC Press.

  EuroMOMO project page, \url{http://www.euromomo.eu/}, Last accessed: 13
  Oct 2010. 
}
\keyword{datasets}
