\name{twinstim_epidataCS_plot}
\alias{plot.epidataCS}

\title{
Plotting the Events of an Epidemic over Time and Space
}

\description{
The \code{plot} method for class \code{"epidataCS"} either plots the
number of events along the time axis (aggregated over space) using
\code{hist()}, or the locations of the events in the observation region
\code{W} (aggregated over time).
}

\usage{
\method{plot}{epidataCS}(x, aggregate = c("time", "space"), subset, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"epidataCS"}.
  }
  \item{aggregate}{
    character, one of \code{"time"} and \code{"space"}. If
    \code{"time"}, the number of events over time is plotted using the
    function \code{\link{hist}} (or \code{\link{hist.Date}}). Otherwise,
    the observation region \code{x$W} and the locations of the events
    therein are plotted.
  }
  \item{subset}{
    logical expression indicating a subset of events to consider for
    plotting: missing values are taken as false. Note that the
    expression is evaluated in the events data frame
    \code{x$events@data}, which means that column names can be referred to
    by name (like in \code{\link{subset.data.frame}}).
  }
  \item{\dots}{
    The following graphical parameters are available depending on the
    value of \code{aggregate}.
    
    For plotting the number of events over time (\code{aggregate = "time"}):
    \describe{
      \item{t0.Date}{the beginning of the observation period
	\code{t0 = x$stgrid$start[1]} as \code{\link{Date}}, enabling
	a nice x-axis using \code{\link{hist.Date}} and sensible
	\code{breaks} of the histogram (which must also be specified in
	this case, e.g., \code{"months"}). The event times then equal
	\code{t0.Date + as.integer(x$events$time - t0)}, i.e. possible
	fractional parts of the event times are removed (which ensures that
	using \code{breaks = "months"} or other automatic types always
	works).}
      \item{freq}{see \code{\link{hist}}, defaults to \code{TRUE}.}
      \item{xlim,ylim}{\code{NULL} provides automatic axis limits.}
      \item{xlab,ylab}{axis labels (with sensible defaults).}
      \item{main}{main title of the plot, defaults to \code{NULL}.}
      \item{col}{fill colour for the bars of the histogram, default to \code{"white"}}
      \item{panel.first}{expression that should be evaluated after the
	plotting window has been set up but before the histogram is plotted.
	Defaults to adding horizontal grid lines.}
      \item{\dots}{further arguments passed to \code{\link{hist}} or
	\code{\link{hist.Date}}, respectively.}
    }

    For plotting the locations of the events in the observation region
    (\code{aggregate = "space"}):
    \describe{
      \item{cex.fun}{function which takes a vector of counts of events
	at each unique location and returns a vector of \code{cex}
	values for the sizes of the corresponding \code{points}.
	Defaults to the \code{sqrt()} function, which for the default
	circular \code{pch=1} means that the area of each point is
	proportional to the number of events at its location.}
      \item{points.args}{a list of further arguments passed to
        \code{\link{points}}. An argument \code{cex} (default: 0.5) in
	this list multiplies the sizes obtained from \code{cex.fun}.}
      \item{colTypes}{colour vector for the different event types.
	Defaults to colours from the \code{\link{rainbow}} palette. They are
	only used if \code{col} is missing in \code{points.args}.}
      \item{\dots}{further arguments passed to the \code{plot}-method
	for \code{"\linkS4class{SpatialPolygons}"}, which draws the
	observation region \code{x$W}.}
    }
  }
}

\value{
For \code{aggregate = "time"} the data of the histogram (as returned by
\code{\link{hist}}), for \code{aggregate = "space"} \code{NULL}, invisibly.
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
plot(imdepi, aggregate = "time", subset = type == "B",
     main = "Time series of events of the serogroup B finetype")

## show the distribution in space
plot(imdepi, aggregate = "space", colTypes = 1:2)
}

\keyword{hplot}
\keyword{methods}
\keyword{spatial}
