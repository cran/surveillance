\name{epidataCS_permute}
\alias{permute.epidataCS}
\title{
  Randomly Permute Time Points or Locations of \code{"epidataCS"}
}

\description{
  Monte Carlo tests for space-time interaction (\code{\link{epitest}})
  use the distribution of some test statistic under the null hypothesis
  of no space-time interaction. For this purpose, the function
  \code{permute.epidataCS} randomly permutes the time or space labels of
  the events.
}

\usage{
permute.epidataCS(x, what = c("time", "space"), keep)
}

\arguments{
  \item{x}{an object of class \code{"\link{epidataCS}"}.}
  \item{what}{character string determining what to permute: time points
    (default) or locations.}
  \item{keep}{optional logical expression to be evaluated in the context
    of \code{x$events@data}, determining for which events the time and
    location should be kept as is. For instance, to keep some
    \dQuote{prehistory} before time point 30 unchanged, use
    \code{keep = time <= 30}.}
}

\value{
  the permuted \code{"\link{epidataCS}"} object.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{epitest}}
}

\examples{
data("imdepi")

set.seed(3)
permepi <- permute.epidataCS(imdepi, what = "time", keep = time <= 30)

print(imdepi, n = 8)
print(permepi, n = 8)
## the first 6 events are kept (as are all row.names),
## the time labels of the remaining events are shuffled
## (and events then again sorted by time),
## the marginal temporal distribution is unchanged
}

\keyword{manip}
