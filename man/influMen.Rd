\name{influMen}
\alias{influMen}
\docType{data}
\title{Influenza and meningococcal infections in Germany, 2001-2006}
\description{
  Weekly counts of new influenza and meningococcal infections in Germany 2001-2006. 
}
\usage{data(influMen)}
\format{
  A \code{disProg} object containing \eqn{312\times 2}{312 x 2}
  observations starting from week 1 in 2001 to week 52 in 2006. 
}
\source{
  Robert Koch-Institut: SurvStat: \url{https://survstat.rki.de/}.
  Queried on 25 July 2007.
}
\examples{
data(influMen)
plot(influMen, as.one=FALSE, same.scale=FALSE) 
}
\keyword{datasets}
