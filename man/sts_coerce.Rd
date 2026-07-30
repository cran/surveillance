\name{sts_conversion}
\alias{as.ts.sts}
\alias{coerce,sts,ts-method}
\alias{coerce,ts,sts-method}
\alias{as.xts.sts}

\title{
  Convert an \code{"sts"} Object to \code{"ts"} or \code{"xts"}
}
\description{
  These methods convert \pkg{surveillance} time series objects of class
  \code{"\linkS4class{sts}"} to the standard \code{"\link{ts}"} and
  \code{"\link[xts]{xts}"} classes.
  Conversion to and from \code{ts} is also possible using S4-type
  \code{coerce} methods, e.g., \code{as(x, "ts")}.
}
\usage{
\method{as.ts}{sts}(x, ...)

\method{as.xts}{sts}(x, order.by = epoch(x, as.Date = TRUE), ...)
}
\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{order.by}{
    optional \code{\link{Date}}s indexing the time series.
    The default is only used for weekly (\code{freq=52}) and daily
    (\code{freq=365}) data. If not specified, other frequencies are
    handled via \code{\link[xts]{as.xts.ts}}.
  }
  \item{\dots}{further arguments are ignored by the \code{as.ts} method
    and passed to \code{\link[xts]{xts}} by the \code{as.xts} method.}
}
\value{
  an object of class \code{"\link{ts}"} or \code{"\link[xts]{xts}"},
  respectively.
}
\examples{
## convert "ts"/"mts" object to "sts"
counts <- matrix(rpois(28*3,10), 28, 3,
                 dimnames = list(NULL, paste("Region", 1:3)))
z <- ts(counts, start = c(1983, 12), frequency = 12)
head(z)
x <- as(z, "sts")
x

## direct "sts" construction gives the same:
(x2 <- sts(counts, start = start(z), frequency = frequency(z)))
\dontshow{stopifnot(identical(x, x2))}

## conversion of "sts" objects to the "xts" class
if (require("xts")) {
    y <- as.xts(x)
    print(head(y))
    plot(y)
}
}

\keyword{ts}
\keyword{methods}
