\encoding{latin1}
\name{twinstim_epidataCS_animate}
\alias{animate.epidataCS}

\title{
Spatio-Temporal Animation of a Continuous-Time Continuous-Space Epidemic
}

\description{
  Function for the animation of continuous-time continuous-space
  epidemic data, i.e. objects inheriting from class \code{"epidataCS"}.
  There are three types of animation, see argument \code{time.spacing}.
  Besides the direct plotting in the \R session, it is 
  also possible to generate a sequence of graphics files to create animations 
  outside \R, since this function is compatible to the framework of the
  contributed \R package \pkg{animation} and respects the
  \code{\link[animation]{ani.options}} \code{"interval"} and
  \code{"nmax"} if this package is available.
}

\usage{
\method{animate}{epidataCS}(object, interval = c(0,Inf), time.spacing = NULL,
        nmax = NULL, sleep = NULL, legend.opts = list(), timer.opts = list(),
        pch = 15:18, col.current = "red", col.I = "#C16E41",
        col.R = "#B3B3B3", col.influence = "#FEE0D2",
        main = NULL, verbose = TRUE, ...)
}

\arguments{
  \item{object}{
    an object inheriting from class \code{"epidataCS"}.
  }
  \item{interval}{time range of the animation.}
  \item{time.spacing}{
    time interval for the animation steps.\cr
    If \code{NULL} (the default), the events are plotted one by one.
    Thus, it is just the \emph{ordering} of the events, which is shown.\cr
    To plot the appearance of events proportionally to the exact time line,
    \code{time.spacing} can be set to a numeric value indicating the period of
    time between consecutive plots.  Then, for each time point in
    \code{seq(0, end, by = time.spacing)} the current state of the epidemic can
    be seen and an additional timer indicates the current time
    (see \code{timer.opts} below).\cr
    If \code{time.spacing = NA}, then the time spacing is automatically
    determined in such a way that \code{nmax} snapshots result. In this
    case, \code{nmax} must be given a finite value.
  }
  \item{nmax}{
    maximum number of snapshots to generate. The default \code{NULL}
    means to take the value from \code{ani.options("nmax")} if the
    \pkg{animation} package is available, and no limitation (\code{Inf})
    otherwise.
  }
  \item{sleep}{
    numeric scalar specyfing the artificial pause in seconds between two
    time points, or \code{NULL} (default), when this is taken from
    \code{ani.options("interval")} if the \pkg{animation} package is
    available, or set to 0.1 otherwise. 
    Thus, if generating an animation inside \pkg{animation}'s
    \code{\link[animation]{saveHTML}} or similar, then this should be
    set to 0.
  }  
  \item{pch, col}{
    vectors of length equal to the number of event types specifying the
    point symbols and colors for events to plot (in this order).
    The vectors are recycled if necessary.
  }
  \item{legend.opts}{
    either a list of arguments passed to the \code{\link{legend}} function or
    \code{NULL} (or \code{NA}), in which case no legend will be plotted.  All 
    necessary arguments have sensible defaults and need not be
    specified.
  }
  \item{timer.opts}{
    either a list of arguments passed to the \code{\link{legend}} function or
    \code{NULL} (or \code{NA}), in which case no timer will be plotted.  All 
    necessary arguments have sensible defaults and need not be specified, i.e.
    \describe{
      \item{\code{x}:}{\code{"bottomright"}}
      \item{\code{title}:}{\code{"time"}}
      \item{\code{box.lty}:}{\code{0}}
      \item{\code{adj}:}{\code{c(0.5,0.5)}}
      \item{\code{inset}:}{\code{0.01}}
      \item{\code{bg}:}{\code{"white"}}
    }
    Note that the argument \code{legend}, which is the current time of the
    animation, can not be modified.
  }
  \item{col.current}{color of events when occuring.}
  \item{col.I}{color once infectious.}
  \item{col.R}{color event has once \dQuote{recovered}. If \code{NA},
    then recovered events will not be shown.}
  \item{col.influence}{color with which the influence region is drawn. Use
    \code{NA} if no influence regions should be drawn.}
  \item{main}{optional main title placed above the map.}
  \item{verbose}{logical specifying if a (textual) progress bar should
    be shown during snapshot generation. This is especially useful if
    the animation is produced within \code{\link[animation]{saveHTML}}
    or similar.}
  \item{\dots}{
    further graphical parameters passed to the plot-method of the
    \code{\link{SpatialPolygons-class}}.
  }
}

%\value{
%  invisibly returns \code{NULL}.
%}

\author{
Sebastian Meyer with documentation contributions by Michael H\enc{ö}{oe}hle
}

\seealso{
\code{\link{plot.epidataCS}} for plotting the numbers of events by time
(aggregated over space) or the locations of the events in the
observation region \code{W} (aggregated over time).

The contributed \R package \pkg{animation}.
}

\examples{
\dontrun{
  data("imdepi")
  imdepiB <- subset(imdepi, type == "B")
  if (require("animation")) ani.options(interval = 0.1, nmax = 60)

  # Animate the first year of type B with a step size of 7 days
  animate(imdepiB, interval=c(0,365), time.spacing=7)

  # Animate the whole time range but with nmax=60 snapshots only and
  # omitting the painting of (quiet artificial) influence regions
  animate(imdepiB, time.spacing=NA, col.influence=NA, nmax=60)

  # Sequential animation of events
  animate(imdepiB, time.spacing=NULL, nmax=60)

  # Such an animation can be saved in various ways using the tools of
  # the animation package, e.g., saveHTML()
  if (require("animation")) {
    ani.options(outdir = tempdir())
    saveHTML(animate(imdepiB, interval=c(0,365), time.spacing=7),
             img.name = "imdepi_year1", htmlfile = "imdepi_animation.html",
             autobrowse = TRUE, title = "Animation of the first year",
             description = c("Take a look into the first-year dynamics."))
  }
}
}

\keyword{hplot}
\keyword{dynamic}
\keyword{spatial}
