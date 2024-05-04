\name{sts_animate}
\alias{animate.sts}
\title{
  Animated Maps and Time Series of Disease Counts or Incidence
}

\description{
  The \code{animate}-method for \code{\linkS4class{sts}} objects
  iterates over time points, plotting maps of the current|cumulative
  counts|incidence via \code{\link{stsplot_space}}, optionally
  including a time series chart below the map to track the epidemic curve.
  It is worth using functionality of the \pkg{animation} package
  (e.g., \code{\link[animation]{saveHTML}}) to directly export the
  animation into a useful format.
}

\usage{
\method{animate}{sts}(object, tps = NULL, cumulative = FALSE,
        population = NULL, at = 10, ...,
        timeplot = list(pos = 1, size = 0.3, fill = TRUE),
        sleep = 0.5, verbose = interactive(), draw = TRUE)
}

\arguments{
  \item{object}{
    an object of class \code{"\linkS4class{sts}"} or a matrix of counts,
    i.e., \code{observed(stsObj)}, where especially
    \code{colnames(x)} have to be contained in \code{row.names(map)}.
    If a matrix, the \code{map} object has to be provided explicitly (as
    part of \code{\dots}).
  }
  \item{tps}{
    a numeric vector of one or more time points at which to plot the map.
    The default \code{tps=NULL} means the whole time period \code{1:nrow(object)}.
  }
  \item{cumulative}{
    logical specifying if the cumulative counts/incidence over time
    should be plotted. The cumulative incidence is relative to the
    population from the first time point \code{tps[1]} throughout the
    whole animation, while \code{cumulative=FALSE} computes the
    incidence from the current population numbers.
  }
  \item{population,at,\dots}{
    arguments for \code{\link{stsplot_space}}.
  }
  \item{timeplot}{
    if a list and package \CRANpkg{gridExtra} is available,
    a time series chart of the counts along
    the selected time points \code{tps} will be plotted next to the map.
    The list elements determine both the positioning of this plot
    (\code{pos}, \code{size}, and \code{fill}) and its appearance.
    The default \code{pos=1} and \code{size=0.3} arguments put the time
    series plot below the map, using 30\% of the total plot height.
    The logical value \code{fill} indicates whether to make the panel as
    big as possible (default: TRUE). An alternative to \code{fill=FALSE}
    is to manually specify an \code{aspect} (ratio) value in \code{timeplot}.
    Other list elements are arguments for the internal (and currently
    undocumented) function \code{stsplot_timeSimple}. For example,
    \code{inactive} and \code{active} are lists of graphical parameters
    (e.g., \code{col}) determining the appearance of the bars (e.g.,
    default color is grey when inactive and black when active),
    and the boolean \code{as.Date} determines whether dates should be
    put on the x-axis (instead of the \code{tps} indexes).
  }
  \item{sleep}{
    time to wait (\code{Sys.sleep}) between subsequent snapshots (only if
    \code{\link{dev.interactive}}), in seconds.
  }
  \item{verbose}{
    logical indicating if a \code{\link{txtProgressBar}} should be shown
    during generation of the animation -- which may take a while.
    Default is to do so in \code{\link{interactive}} sessions.
  }
  \item{draw}{
    logical indicating if the produced plots at each time point should
    be drawn directly (the default) or not.
    The setting \code{draw = FALSE} is useful if one would like to
    manually arrange the plots, which are always returned invisibly in a
    list of length \code{length(tps)}.
  }
}

\value{
  (invisibly) a list of the \code{length(tps)} sequential plot objects.
  These are of class \code{"gtable"} (from \CRANpkg{gtable})
  if the \code{timeplot} is included, otherwise of class
  \code{"\code{\link[lattice:trellis.object]{trellis}"}.
}

\author{
  Sebastian Meyer
}

\seealso{
  the other plot types documented in \code{\link{stsplot}} for static
  time series plots and maps.
}

\examples{
data("measlesWeserEms")

## animate the weekly counts of measles (during weeks 12-16 only, for speed)
if (interactive() && require("animation")) {
    oldwd <- setwd(tempdir())  # to not clutter up the current working dir
    saveHTML(animate(measlesWeserEms, tps=12:16),
             title="Evolution of the measles epidemic in the Weser-Ems region",
             ani.width=500, ani.height=600)
    setwd(oldwd)
}

\dontrun{
## animate the weekly incidence of measles (per 100'000 inhabitants),
## and label the time series plot with dates in a specified format
animate(measlesWeserEms, tps=12:16,
        population = measlesWeserEms@map$POPULATION / 100000,
        timeplot = list(as.Date = TRUE,
                        scales = list(x = list(format = "\%G/\%V"))))
}
}

\keyword{hplot}
\keyword{dynamic}
\keyword{spatial}
