\name{formatDate}
\alias{formatDate}
\title{
Convert Dates to Character (Including Quarter Strings)
}
\description{
  An extension of \code{\link{format.Date}} with additional formatting
  strings for quarters. Used by \code{\link{linelist2sts}}.
}
\usage{
formatDate(x, format)
}
\arguments{
  \item{x}{a \code{"\link{Date}"} object.}
  \item{format}{
    a character string, see \code{\link{strftime}} for possible
    specifications. Further to these base formats, \code{formatDate}
    implements:
    \describe{
      \item{\code{"\%Q"}}{the quarter as a numeric}
      \item{\code{"\%OQ"}}{the quarter as a roman numeral}
      \item{\code{"\%q"}}{the day within the quarter}
    }
  }
}
\value{
  a character vector representing the input date(s) \code{x}
  following the \code{format} specification.
}
\seealso{
\code{\link{strftime}}
}

\examples{
formatDate(Sys.Date(), "\%G/\%OQ/\%q")
}

\keyword{chron}
