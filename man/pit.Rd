\name{pit}
\alias{pit}

\title{
Non-Randomized Version of the PIT Histogram (for Count Data)
}

\description{
  See Czado et al. (2009).
}

\usage{
pit(x, pdistr, J = 10, relative = TRUE, ..., plot = list())
}

\arguments{
  \item{x}{
    numeric vector representing the observed counts.
  }
  \item{pdistr}{
    either a list of predictive cumulative distribution functions for
    the observations \code{x}, or (the name of) a single predictive CDF
    used for all \code{x} (with potentially varying arguments \code{...}).
    It is checked that the predictive CDF returns 0 at \code{x=-1}.
    The name of its first argument can be different from \code{x}, e.g.,
    \code{pdistr="pnbinom"} is possible.\cr
    If \code{pdistr} is a single function and no additional \code{\dots}
    arguments are supplied, \code{pdistr} is assumed to be vectorized,
    i.e., it is simply called as \code{pdistr(x)} and \code{pdistr(x-1)}.
    Otherwise, the predictive CDF is called sequentially and does not
    need to be vectorized.
  }
  \item{J}{
    the number of bins of the histogram.
  }
  \item{relative}{
    logical indicating if relative frequency or the density should be plotted.
  }
  \item{\dots}{
    ignored if \code{pdistr} is a list. Otherwise, such additional
    arguments are used in sequential calls of \code{pdistr} via
    \code{\link{mapply}(pdistr, x, ...)}.
  }
  \item{plot}{
    a list of arguments for \code{\link{plot.histogram}} or \code{NULL}
    in which case no plot will be produced.
  }
}

\value{
  an object of class \code{"histogram"} (see \code{\link{hist}}).
  It is returned invisibly if a plot is produced.
}

\references{
  Czado, C., Gneiting, T. & Held, L. (2009):
  Predictive model assessment for count data.
  \emph{Biometrics}, \bold{65}, 1254-1261.
}

\author{
  Michaela Paul and Sebastian Meyer
}

\examples{
## Simulation example of Czado et al. (2009, Section 2.4)
set.seed(100)
x <- rnbinom(200, mu = 5, size = 2)
pdistrs <- list("NB(5,0)"   = function (x) ppois(x, lambda=5),
                "NB(5,1/2)" = function (x) pnbinom(x, mu=5, size=2),
                "NB(5,1)"   = function (x) pnbinom(x, mu=5, size=1))
## Reproduce Figure 1
op <- par(mfrow = c(1,3))
for (i in seq_along(pdistrs)) {
    pit(x, pdistr = pdistrs[[i]], J = 10, relative = TRUE,
        plot = list(ylim = c(0,2.75), main = names(pdistrs)[i]))
    box()
}
par(op)

## Alternative call using ... arguments for pdistr (less efficient)
stopifnot(identical(pit(x, "pnbinom", mu = 5, size = 2, plot = NULL),
                    pit(x, pdistrs[[2]], plot = NULL)))
}

\keyword{dplot}