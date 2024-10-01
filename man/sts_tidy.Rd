\name{tidy.sts}
\alias{tidy.sts}
\title{
  Convert an \code{"sts"} Object to a Data Frame in Long (Tidy) Format
}

\description{
  The resulting data frame will have a row for each time point and
  observational unit, and columns corresponding to the slots of the
  \code{"\linkS4class{sts}"} object (except for \code{populationFrac},
  which is named \code{population}).
  Some time variables are added for convenience:
  \code{year}, \code{epochInYear}, \code{epochInPeriod}, \code{date}
  (the latter gives \code{NA} dates if \code{epoch(x, as.Date=TRUE)}
  fails, i.e., for non-standard \code{frequency(x)} if \code{x@epochAsDate} is false).
}

\usage{
tidy.sts(x, ...)
}

\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{\dots}{unused.}
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{as.data.frame.sts}}
}

\examples{
data("momo")
momodat <- tidy.sts(momo)
head(momodat)

## tidy.sts(stsObj) is the same as as.data.frame(stsObj, tidy = TRUE)
stopifnot(identical(as.data.frame(momo, tidy = TRUE), momodat))
}

\keyword{manip}
