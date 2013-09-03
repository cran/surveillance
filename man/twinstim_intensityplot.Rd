\name{twinstim_intensityplot}
\alias{intensityplot.twinstim}
\alias{intensity.twinstim}
\alias{intensityplot.simEpidataCS}
\title{
  Plotting Intensities of Infection over Time or Space
}

\description{
\code{\link{intensityplot}} method to plot the evolution of the total infection
intensity, its epidemic proportion or its endemic proportion over time
or space (integrated over the other dimension) of fitted
\code{\link{twinstim}} models (or \code{\link{simEpidataCS}}).
The \code{"simEpidataCS"}-method is just a wrapper around
\code{intensityplot.twinstim} by making the \code{"simEpidataCS"} object
\code{"twinstim"}-compatible, i.e. enriching it by the
required model components and environment.

The \code{intensity.twinstim} auxiliary function returns functions which
calculate the endemic or epidemic intensity at a specific time point or
location (integrated over the other dimension).
}

\usage{
\method{intensityplot}{twinstim}(x,
    which = c("epidemic proportion", "endemic proportion", "total intensity"),
    aggregate = c("time", "space"), types = 1:nrow(x$qmatrix),
    tiles, tiles.idcol = NULL, plot = TRUE, add = FALSE,
    tgrid = 101, rug.opts = list(),
    sgrid = 128, polygons.args = list(), points.args = list(),
    cex.fun = sqrt, ...)

\method{intensityplot}{simEpidataCS}(x, ...)

intensity.twinstim(x,
    aggregate = c("time", "space"), types = 1:nrow(x$qmatrix), 
    tiles, tiles.idcol = NULL)
}

\arguments{
  \item{x}{
    an object of class \code{"twinstim"} or \code{"simEpidataCS"}, respectively.
}
  \item{which}{
    \code{"epidemic proportion"}, \code{"endemic proportion"},
    or \code{"total intensity"}.  Partial matching is applied.  Determines
    whether to plot the path of the total intensity or its 
    epidemic or endemic proportions over time or space (\code{which})
    aggregated over the other dimension and \code{types}.
}
  \item{aggregate}{
    One of \code{"time"} or \code{"space"}. The former results in a plot
    of the evolution of \code{which} as a function of time (integrated
    over the observation region \eqn{\bold{W}}), whereas the latter
    produces a \code{spplot} of \code{which} over \eqn{\bold{W}}
    (spanned by \code{tiles}). In both cases, \code{which} is evaluated
    on a grid of values, given by \code{tgrid} or \code{sgrid}, respectively.
}
  \item{types}{
    event types to aggregate. By default, all types of events are
    aggregated, but one could also be interested in only one specific
    type or a subset of event types.
}
  \item{tiles}{
    object of class \code{\linkS4class{SpatialPolygons}} representing
    the decomposition of \eqn{\bold{W}} into different regions (as used
    in the corresponding \code{stgrid} of the \code{"\link{epidataCS}"}.
    This is only needed for \code{aggregate = "space"}.
}
  \item{tiles.idcol}{
    either a column index for \code{tiles@data}, or \code{NULL}
    (default), which refers to the \code{"ID"} slot of the polygons.
    The resulting ID's must correspond to the ID's used in the
    \code{tile} column of the \code{stgrid} of the
    \code{"\link{epidataCS}"} on which \code{x} was fitted.
}
  \item{plot}{
    logical indicating if a plot is desired, which defaults to \code{TRUE}.
    Otherwise, a function will be returned, which takes a vector of time
    points (if \code{aggregate = "time"}) or a matrix of coordinates (if
    \code{aggregate = "space"}), and returns \code{which} on this grid.
  }
  \item{add}{
    logical.  If \code{TRUE} and \code{aggregate = "time"}, paths are
    added to the current plot, using \code{lines}.  This does not work
    for \code{aggregate = "space"}.
  }
  \item{tgrid}{
    either a numeric vector of time points when to evaluate
    \code{which}, or a scalar representing the desired number of
    evaluation points in the observation interval \eqn{[t_0, T]}.
    This argument is unused for \code{aggregate = "space"}.
  }
  \item{rug.opts}{
    if a list, its elements are passed as arguments to the function
    \code{\link{rug}}, which will mark the time points of the events if
    \code{aggregate = "time"} (it is unused in the spatial case);
    otherwise (e.g., \code{NULL}), no \code{rug} will be produced.
    By default, the \code{rug} argument \code{ticksize} is set to 0.02
    and \code{quiet} is set to \code{TRUE}.  Note that the argument
    \code{x} of the \code{rug} function, which contains the
    locations for the \code{rug} is fixed internally and can not be
    modified.
  }
  \item{sgrid}{
    either an object of class \code{"\linkS4class{SpatialPixels}"} (or
    coercible to that class) representing the locations where to
    evaluate \code{which}, or a scalar representing the total
    number of points of a grid constructed on the bounding box of
    \code{tiles} (which requires the package \pkg{maptools}).
    \code{sgrid} is internally subsetted to contain only
    points inside \code{tiles}.
    This argument is unused for \code{aggregate = "time"}.
}
  \item{polygons.args}{
    if a list, its elements are passed as arguments to
    \code{\link{sp.polygons}}, which will add \code{tiles} to the plot
    if \code{aggregate = "space"} (it is unused for the temporal plot).
    By default, the fill \code{col}our of the tiles is set to
    \code{"darkgrey"}.
}
  \item{points.args}{
    if a list, its elements are passed as arguments to
    \code{\link{sp.points}}, which will add the event locations to the plot
    if \code{aggregate = "space"} (it is unused for the temporal plot).
    By default, the plot symbol ist set to \code{pch=1}. The sizes
    of the points are determined as the product of the argument \code{cex}
    (default: 0.5) of this list and the sizes obtained from
    the function \code{cex.fun} which accounts for multiple events at the
    same location.
}
  \item{cex.fun}{
    function which takes a vector of counts of events
    at each unique location and returns a (vector of) \code{cex}
    value(s) for the sizes of the points at the event locations used in
    \code{points.args}.
    Defaults to the \code{sqrt()} function, which for the default
    circular \code{pch=1} means that the area of each point is
    proportional to the number of events at its location.
}
  \item{\dots}{
    further arguments passed to \code{plot} or \code{lines} (if
    \code{aggregate = "time"}), or to \code{\link{spplot}} (if
    \code{aggregate = "space"}).\cr
    For \code{intensityplot.simEpidataCS}, arguments passed to
    \code{intensityplot.twinstim}.
}
}

\value{
  If \code{plot = FALSE} or \code{aggregate = "time"},
  a function is returned, which takes a vector of
  time points (if \code{aggregate = "time"}) or a matrix of coordinates
  (if \code{aggregate = "space"}), and returns \code{which} on this grid.

  If \code{plot = TRUE} and \code{aggregate = "space"}, the
  \code{\link{trellis.object}} containg the spatial plot is returned.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{plot.twinstim}}, which calls this function.
}

\examples{
data("imdepi")
data("imdepifit")

# for the intensityplot we need the model environment, which can be
# easily added by the intelligent update method (no need to refit the model)
imdepifit <- update(imdepifit, model=TRUE)

## path of the total intensity
opar <- par(mfrow=c(2,1))
intensityplot(imdepifit, which="total intensity",
              aggregate="time", tgrid=500)
plot(imdepi, breaks=100)
par(opar)

## path of the proportion of the epidemic component by event
intensityplot(imdepifit, which="epidemic proportion",
              aggregate="time", tgrid=500, types=1)
intensityplot(imdepifit, which="epidemic proportion",
              aggregate="time", tgrid=500, types=2, add=TRUE, col=2)
legend("topright", legend=levels(imdepi$events$type), lty=1, col=1:2,
       title = "event type")

if (surveillance.options("allExamples") && requireNamespace("maptools"))
{
  ## spatial shape of the intensity (aggregated over time)
  # load borders of Germany's districts (obtained from the
  # Bundesamt für Kartographie und Geodäsie, Frankfurt am Main, Germany,
  # www.geodatenzentrum.de), simplified by the "special" Visvalingam
  # algorithm (level=60\%) using www.MapShaper.org:
  load(system.file("shapes", "districtsD.RData", package="surveillance"))

  # total intensity (using a rather sparse 'sgrid' for speed)
  intensityplot(imdepifit, which="total intensity",
                aggregate="space", tiles=districtsD, sgrid=500)

  # epidemic proportion by type
  maps_epiprop <- lapply(1:2, function (type) {
      intensityplot(imdepifit, which="epidemic", aggregate="space",
                    types=type, tiles=districtsD, sgrid=1000,
                    at=seq(0,1,by=0.1), col.regions=rev(heat.colors(20)))
  })
  plot(maps_epiprop[[1]], split=c(1,1,2,1), more=TRUE)
  plot(maps_epiprop[[2]], split=c(2,1,2,1))
}
}

\keyword{hplot}
\keyword{aplot}
\keyword{dplot}
\keyword{methods}
