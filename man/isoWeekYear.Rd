\name{isoWeekYear}
\alias{isoWeekYear}

\title{Find ISO Week and Year of Date Objects}
\description{
  The function \code{isoWeekYear} extracts the year and week of a
  \code{\link{Date}} according to the ISO 8601 specification. It simply
  calls \code{\link{strftime}} with format strings \code{"\%G"}
  and \code{"\%V"}, respectively.
}
\usage{
isoWeekYear(Y, M, D)
}
\arguments{
  \item{Y}{year(s) or a Date/POSIXt object. Can be a vector.}
  \item{M}{month(s), only used if \code{Y} is not a Date/POSIXt object.}
  \item{D}{day(s), only used if \code{Y} is not a Date/POSIXt object.}
}
\value{
  A list with entries \code{ISOYear} and \code{ISOWeek} containing the
  corresponding results.
}
\examples{
dates <- as.Date(c("2002-12-31","2003-01-01","2003-01-06"))
isoWeekYear(dates)

## the same using numeric inputs:
isoWeekYear(Y = c(2002, 2003, 2003), M = c(12, 1, 1), D = c(31, 1, 6))
}

\keyword{chron}
