\name{sts_ggplot}
\alias{autoplot.sts}
\title{
  Time-Series Plots for \code{"sts"} Objects Using \pkg{ggplot2}
}
\description{
  A simple \CRANpkg{ggplot2} variant of \code{\link{stsplot_time}},
  based on a \dQuote{tidy} version of the \code{"sts"} object via
  \code{\link{tidy.sts}}.
}
\usage{
autoplot.sts(object, population = FALSE, units = NULL,
             as.one = FALSE, scales = "fixed", width = NULL, ...)
}
\arguments{
  \item{object}{an object of class \code{"\linkS4class{sts}"}.}
  \item{population}{logical indicating whether \code{observed(object)}
    should be divided by \code{population(object)}.
    The \code{population} argument can also be a scalar,
    which is used to scale the denominator \code{population(object)},
    i.e., \code{observed(object)} is divided by
    \code{population(object) / population}.
    For instance, if \code{population(object)} contains raw population
    numbers, \code{population = 1000} could be used to plot the
    incidence per 1000 inhabitants.}
  \item{units}{optional integer or character vector to select the units
    (=columns of \code{object}) to plot. The default (\code{NULL}) is
    to plot all time series.}
  \item{as.one}{logical indicating if all time series should be plotted
    in one panel with \code{\link[ggplot2]{geom_line}}.
    By default, the time series are plotted in separate panels (using
    \code{\link[ggplot2]{geom_col}}).}
  \item{scales}{passed to \code{\link[ggplot2]{facet_wrap}}
    (for \code{as.one=FALSE}). By default, all panels use a common
    \code{ylim} (and \code{xlim}).}
  \item{width}{bar width, passed to \code{\link[ggplot2]{geom_col}}.
    Defaults to 7 for weekly time series.}
  \item{\dots}{unused (argument of the generic).}
}
\value{
  a \code{"ggplot"} object.
}
\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{stsplot_time}} for the traditional plots.
}

\examples{
## compare traditional plot() with ggplot2-based autoplot.sts()
if (requireNamespace("ggplot2")) {
    data("measlesDE")
    plot(measlesDE)
    autoplot.sts(measlesDE)
}

## weekly incidence: population(measlesDE) gives population fractions,
## which we need to multiply by the total population
if (surveillance.options("allExamples") && require("ggplot2")) {
    autoplot.sts(measlesDE, population = 1000000/82314906) +
        ylab("Weekly incidence [per 1'000'000 inhabitants]")
}
}

\keyword{hplot}
\keyword{ts}
