\name{salmHospitalized}
\alias{salmHospitalized}
\docType{data}
\title{Hospitalized Salmonella cases in Germany 2004-2014}
\description{
  Reported number of cases of Salmonella in Germany 2004-2014 (early 2014) that were hospitalized. The corresponding 
  total number of cases is indicated in the \code{slot populationFrac}. \code{multinomialTS} is \code{TRUE}.
}
\usage{data(salmNewport)}
\format{
  A \code{sts} object. 
}
\source{  The data are queried from the Survstat@RKI database of the German
  Robert Koch Institute (http://www3.rki.de/SurvStat/). }
\examples{
data(salmHospitalized)
}
\keyword{datasets}
