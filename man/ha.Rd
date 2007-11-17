\name{ha}
\alias{ha}
\title{Hepatitis A in Berlin}
\description{Number of Hepatitis A cases among adult male (age>18) in
  Berlin 2001-2006. An increase is seen during 2006
}
\usage{data(ha)
}
\format{
  A \code{disProg} object containing \eqn{290\times 12}{290 x 12}
  observations starting from week 1 in 2001 to week 30 in 2006. 
}
\source{
Robert Koch-Institut: SurvStat: \url{http://www3.rki.de/SurvStat};
Queried on 25 August 2006.

Robert Koch Institut, Epidemiologisches Bulletin 33/2006, p.290.
}
\examples{
data(ha)
plot(aggregate(ha))
}
\keyword{datasets}


