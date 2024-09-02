\name{shadar}
\alias{shadar}
\docType{data}
\title{Salmonella Hadar cases in Germany 2001-2006}
\description{
  Number of salmonella hadar cases in
  Germany 2001-2006. An increase is seen during 2006.
}
\usage{data(shadar)}
\format{
  A \code{disProg} object containing \eqn{295\times 1}{295 x 1}
  observations starting from week 1 in 2001 to week 35 in 2006. 
}
\source{
Robert Koch-Institut: SurvStat: \url{https://survstat.rki.de/};
Queried on September 2006.

Robert Koch Institut, Epidemiologisches Bulletin 31/2006.
}
\examples{
data(shadar)
plot(shadar)
}
\keyword{datasets}
