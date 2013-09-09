\name{intensityplot}
\alias{intensityplot}
\title{
  Plot Paths of Point Process Intensities
}

\description{
  Generic function for plotting paths of point process intensities.
  Methods currently defined in package \pkg{surveillance} are for
  classes \code{"twinSIR"} and \code{"simEpidata"} (temporal), as well as
  \code{"twinstim"} and \code{"simEpidataCS"} (spatio-temporal). 
}

\usage{
intensityplot(x, ...)
}

\arguments{
  \item{x}{
    An object for which an \code{intensityplot} method is defined.
  }
  \item{\dots}{
    Arguments passed to the corresponding method.
  }
}

\seealso{
  The methods \code{\link{intensityplot.twinSIR}} and
  \code{\link{intensityplot.twinstim}}.
}

\keyword{hplot}
