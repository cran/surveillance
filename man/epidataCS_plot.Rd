\name{epidataCS_plot}
\alias{plot.epidataCS}
\alias{epidataCSplot_time}
\alias{epidataCSplot_space}

\title{
Plotting the Events of an Epidemic over Time and Space
}

\description{
The \code{plot} method for class \code{"epidataCS"} either plots the
number of events along the time axis (\code{epidataCSplot_time}) as a
\code{hist()}, or the locations of the events in the observation region
\code{W} (\code{epidataCSplot_space}).
}

\usage{
\method{plot}{epidataCS}(x, aggregate = c("time", "space"), subset, ...)

epidataCSplot_time(x, subset, t0.Date = NULL, freq = TRUE,
                   col = rainbow(nTypes), cumulative = list(),
                   add = FALSE, mar = NULL, xlim = NULL, ylim = NULL,
                   xlab = "Time", ylab = NULL, main = NULL,
                   panel.first = abline(h=axTicks(2), lty=2, col="grey"),
                   legend.types = list(), ...)

epidataCSplot_space(x, subset,
                    cex.fun = sqrt, points.args = list(), add = FALSE,
                    legend.types = list(), legend.counts = list(), ...)
}

\arguments{
  \item{x}{
    an object of class \code{"\link{epidataCS}"}.
  }
  \item{aggregate}{
    character, one of \code{"time"} and \code{"space"}, referring to the
    specific plot functions \code{epidataCSplot_time} and
    \code{epidataCSplot_time}, respectively.
    For \code{"time"}, the number of events over time is plotted as
    \code{\link{hist}} (or \code{\link{hist.Date}}).
    For \code{"space"}, the observation region \code{x$W} and the
    locations of the events therein are plotted.
  }
  \item{subset}{
    logical expression indicating a subset of events to consider for
    plotting: missing values are taken as false. Note that the
    expression is evaluated in the data frame of event marks
    (\code{marks(x)}), which means that column names can be referred to
    by name (like in \code{\link{subset.data.frame}}).
  }
  \item{\dots}{
    in the basic \code{plot}-method further arguments are passed to the
    \code{aggregate}-specific plot function.
    In \code{epidataCSplot_time}, further graphical parameters are
    passed to \code{\link{hist}} or \code{\link{hist.Date}},
    respectively. In \code{epidataCSplot_space}, further arguments are
    passed to the \code{plot}-method for
    \code{"\linkS4class{SpatialPolygons}"}, which draws the observation
    region \code{x$W}.
  }
    
  \item{t0.Date}{the beginning of the observation period
    \code{t0 = x$stgrid$start[1]} as a \code{"\link{Date}"} (or
    anything coercible by \code{as.Date} without further arguments),
    enabling a nice x-axis using \code{\link{hist.Date}} and sensible
    \code{breaks} of the histogram (which must also be specified in
    this case, e.g., \code{breaks="months"}). The event times then equal
    \code{t0.Date + as.integer(x$events$time - t0)}, i.e. possible
    fractional parts of the event times are removed (which ensures that
    using \code{breaks = "months"} or other automatic types always
    works).}
  \item{freq}{see \code{\link{hist}}, defaults to \code{TRUE}.}
  \item{col}{fill colour for the bars of the histogram, defaults to
    the vector of \code{\link{rainbow}} colours.}
  \item{cumulative}{if a list (of style options),
    lines for the cumulative number of events (per type) will be
    added to the plot. Possible options are \code{axis} (logical),
    \code{lab} (axis label), \code{maxat} (single integer affecting
    the axis range), \code{lwd}, \code{col}, and \code{offset} (a
    numeric vector of length the number of types).}
  \item{add}{logical (default: \code{FALSE}) indicating if the plot
    should be added to an existing window.}
  \item{mar}{see \code{\link{par}}. The default (\code{NULL}) is
    \code{mar <- par("mar")}, with \code{mar[4] <- mar[2]} if an
    axis is requested for the \code{cumulative} numbers.}
  \item{xlim,ylim}{\code{NULL} provides automatic axis limits.}
  \item{xlab,ylab}{axis labels (with sensible defaults).}
  \item{main}{main title of the plot (defaults to no title).}
  \item{panel.first}{expression that should be evaluated after the
    plotting window has been set up but before the histogram is plotted.
    Defaults to adding horizontal grid lines.}
  \item{legend.types}{if a list (of arguments for \code{\link{legend}}),
    a legend for the event types is added to the plot.}

  \item{cex.fun}{function which takes a vector of counts of events
    at each unique location and returns a (vector of) \code{cex}
    value(s) for the sizes of the corresponding \code{points}.
    Defaults to the \code{sqrt()} function, which for the default
    circular \code{pch=1} means that the area of each point is
    proportional to the number of events at its location.}
  \item{points.args}{a list of (type-specific) graphical parameters
    for \code{\link{points}}, specifically \code{pch}, \code{lwd},
    and \code{col}, which are all recycled to give the length
    \code{nlevels(x$events$type)}. In contrast, a possible
    \code{cex} element should be scalar (default: 0.5) and 
    multiplies the sizes obtained from \code{cex.fun}.}
  \item{legend.counts}{if a list (of arguments for
    \code{\link{legend}}), a legend illustrating the effect of 
    \code{cex.fun} is added to the plot. This list may contain a
    special element \code{counts}, which is an integer vector
    specifying the counts to illustrate.}
}

\value{
For \code{aggregate="time"} (i.e., \code{epidataCSplot_time}) the data
of the histogram (as returned by \code{\link{hist}}),
and for \code{aggregate="space"} (i.e., \code{epidataCSplot_space})
\code{NULL}, invisibly.
}

\author{
Sebastian Meyer
}

\seealso{
\code{\link{animate.epidataCS}}
}

\examples{
data("imdepi")

## show the occurrence of events along the time axis (-> histogram)
plot(imdepi, "time", main = "Histogram of event time points")

## show the distribution in space
plot(imdepi, "space", lwd=2)
}

\keyword{hplot}
\keyword{methods}
\keyword{spatial}
