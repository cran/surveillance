\name{twinstim_plot}
\alias{plot.twinstim}

\title{
Plot methods for fitted \code{twinstim}'s
}

\description{
The fitted conditional intensity function from \code{\link{twinstim}}
may be visualized in at least two ways: \code{\link{iafplot}} plots the fitted
interaction functions (as a function of the distance from the host), and
\code{\link{intensityplot.twinstim}} plots the fitted intensity either
aggregated over space (evolution over time) or aggregated over time
(spatial surface of the cumulated intensity). The \code{plot} method for
class \code{"twinstim"} is just a wrapper for these two functions.
}

\usage{
\method{plot}{twinstim}(x, which, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"twinstim"}.
  }
  \item{which}{
    character. Which characteristic of the conditional intensity should
    be plotted? Possible values are the ones allowed in
    the functions \code{\link{iafplot}} and
    \code{\link{intensityplot.twinstim}}, e.g. \code{"siaf"}, or
    \code{"epidemic proportion"}. Partial matching is applied.
  }
  \item{\dots}{
    further arguments passed to \code{iafplot} or
    \code{intensityplot.twinstim}.
  }
}

\value{
See the documentation of the respective plot functions,
\code{\link{iafplot}} or \code{\link{intensityplot.twinstim}}.
}

\author{
Sebastian Meyer
}

\examples{
# see the examples for iafplot() and intensityplot.twinstim()
}

\keyword{hplot}
