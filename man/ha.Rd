\name{ha}
\alias{ha}
\alias{ha.sts}
\docType{data}
\title{Hepatitis A in Berlin}

\description{
  Number of Hepatitis A cases among adult (age>18) males in
  Berlin, 2001-2006. An increase is seen during 2006.
}

\usage{
data("ha")
data("ha.sts")
}

\format{
  \code{ha} is a \code{disProg} object containing \eqn{290\times 12}{290 x 12}
  observations starting from week 1 in 2001 to week 30 in 2006.
  \code{ha.sts} was generated from \code{ha} via the converter function
  \code{\link{disProg2sts}} and includes a map of Berlin's districts.
}

\source{
Robert Koch-Institut: SurvStat: \url{https://survstat.rki.de/};
Queried on 25 August 2006.

Robert Koch Institut, Epidemiologisches Bulletin 33/2006, p.290.
}

\examples{
## deprecated "disProg" object
data("ha")
ha
plot(aggregate(ha))

## new-style "sts" object
data("ha.sts")
ha.sts
plot(ha.sts, ~time)  # = plot(aggregate(ha.sts, by = "unit"))
plot(ha.sts, ~unit, labels = TRUE)
}

\keyword{datasets}
