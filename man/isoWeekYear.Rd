\name{isoWeekYear}
\alias{isoWeekYear}
\alias{formatDate}

\title{Find ISO week and ISO year of a vector of Date objects on Windows}
\description{
  This function extracts the ISO week and ISO year of a \code{Date}
  according to the ISO 8601 specification. Note that this function
  does nothing else than format.Date(x, "\%G") and format.Date(x, "\%V")
  would do on Mac/Unix computers. However, this is not implemented
  on Windows.

  A small internal wrapper for \code{format.Date} (called
  \code{formatDate}) thus directs all calls having one of these format
  strings to this function, if the \code{sessionInfo()[[1]]$os}
  information reveals a Windows system. The function also provides three
  additional \code{strptime} formatting strings: "\%Q" (the quarter of a
  date as a numeric), "%OQ" (the quarter of a date in roman) and "%q"
  (day within the quarter). These are, e.g., used by \code{linelist2sts}.
  
}
\details{
  The code to find the ISO week and year on Windows is by Gustaf Rydevik
  <gustaf.rydevik_at_gmail.com>  posted
  at http://tolstoy.newcastle.edu.au/R/e10/help/10/05/5588.html
}
\usage{
isoWeekYear(Y, M=NULL, D=NULL)
}
\arguments{
  \item{Y}{Date object (POSIX) or the year. Can be a vector.}
  \item{M}{month, \code{NULL} if Y is a Date object)}
  \item{D}{day, \code{NULL} if Y is a Date object)}
}
\value{
  A list with entries \code{ISOYear} and \code{ISOWeek} containing the
  corresponding results.
}
\author{Gustaf Rydevik}
\examples{
dates <- as.Date(c("2002-12-31","2003-01-01","2003-01-06"))
isoWeekYear(dates)
}

\keyword{chron}

