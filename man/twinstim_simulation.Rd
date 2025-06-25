\encoding{latin1}
\name{twinstim_simulation}
\alias{simEpidataCS}
\alias{simulate.twinstim}

\title{
  Simulation of a Self-Exciting Spatio-Temporal Point Process
}

\description{
  The function \code{simEpidataCS} simulates events of a self-exciting
  spatio-temporal point process of the \code{"\link{twinstim}"} class.
  Simulation works via Ogata's modified thinning of the conditional
  intensity as described in Meyer et al. (2012). Note that simulation is
  limited to the spatial and temporal range of \code{stgrid}.

  The \code{\link{simulate}} method for objects of class
  \code{"\link{twinstim}"} simulates new epidemic data using the model and
  the parameter estimates of the fitted object.
}

\usage{
simEpidataCS(endemic, epidemic, siaf, tiaf, qmatrix, rmarks,
    events, stgrid, tiles, beta0, beta, gamma, siafpars, tiafpars,
    epilink = "log", t0 = stgrid$start[1], T = tail(stgrid$stop,1),
    nEvents = 1e5, control.siaf = list(F=list(), Deriv=list()),
    W = NULL, trace = 5, nCircle2Poly = 32, gmax = NULL, .allocate = 500,
    .skipChecks = FALSE, .onlyEvents = FALSE)

\method{simulate}{twinstim}(object, nsim = 1, seed = NULL, data, tiles,
    newcoef = NULL, rmarks = NULL, t0 = NULL, T = NULL, nEvents = 1e5,
    control.siaf = object$control.siaf,
    W = data$W, trace = FALSE, nCircle2Poly = NULL, gmax = NULL,
    .allocate = 500, simplify = TRUE, ...)
}

\arguments{
  \item{endemic}{
    see \code{\link{twinstim}}. Note that type-specific endemic
    intercepts are specified by \code{beta0} here, not by the term
    \code{(1|type)}.
  }
  \item{epidemic}{
    see \code{\link{twinstim}}. Marks appearing in this formula must
    be returned by the generating function \code{rmarks}.
  }
  \item{siaf}{
    see \code{\link{twinstim}}.
    In addition to what is required for fitting with \code{twinstim},
    the \code{siaf} specification must also contain the element
    \code{simulate}, a function which draws random locations following the
    spatial kernel \code{siaf$f}. The first argument of the function is the
    number of points to sample (say \code{n}),
    the second one is the vector of parameters
    \code{siafpars}, the third one is the type indicator (a character string
    matching a type name as specified by \code{dimnames(qmatrix)}). With the
    current implementation there will always be simulated only one
    location at a time, i.e. \code{n=1}.
    The \link[=siaf.constant]{predefined siaf's} all provide simulation.
  } 
  \item{tiaf}{
    e.g. what is returned by the generating function
    \code{\link{tiaf.constant}} or \code{\link{tiaf.exponential}}. See also
    \code{\link{twinstim}}.
  }
  \item{qmatrix}{
    see \code{\link{epidataCS}}. Note that this square
    matrix and its \code{dimnames} determine the number and names of the
    different event types. In the simplest case, there is only a single
    type of event, i.e. \code{qmatrix = diag(1)}.
  }
  \item{rmarks}{
    function of single time (1st argument) and location
    (2nd argument) returning a one-row \code{data.frame} of marks (named
    according to the variables in \code{epidemic}) for an event at this
    point. This must include the columns \code{eps.s} and \code{eps.t}, 
    i.e. the values of the spatial and temporal interaction ranges at this
    point. Only \code{"numeric"} and \code{"factor"} columns are
    allowed. Assure that factor variables are coded equally 
    (same levels and level order) for each new sample.

    For the \code{simulate.twinstim} method, the default (\code{NULL})
    means sampling from the empirical distribution function of the
    (non-missing) marks in \code{data} restricted to events in the
    simulation period (\code{t0};\code{T}]. If there are no events in
    this period, e.g., if simulating beyond the original observation
    period, \code{rmarks} will sample marks from all of
    \code{data$events}.
  }
  \item{events}{
    \code{NULL} or missing (default) in case of an empty prehistory,
    or a \code{\link[sp]{SpatialPointsDataFrame}} containing events of the
    prehistory (-Inf;\code{t0}] of the process (required for the
    epidemic to start in case of no endemic component in the model).
    The \code{SpatialPointsDataFrame} must have the same
    \code{proj4string} as \code{tiles} and \code{W}). The attached
    \code{data.frame} (data slot) must contain the typical 
    columns as described in \code{\link{as.epidataCS}} (\code{time},
    \code{tile}, \code{eps.t}, \code{eps.s}, and, for type-specific
    models, \code{type}) and all marks appearing in the \code{epidemic}
    specification. Note that some column names are reserved (see
    \code{\link{as.epidataCS}}).  Only events up to
    time \code{t0} are selected and taken as the prehistory.
  }
  \item{stgrid}{
    see \code{\link{as.epidataCS}}. Simulation only works inside the spatial
    and temporal range of \code{stgrid}.
  }
  \item{tiles}{
    object inheriting from \code{"\linkSPclass{SpatialPolygons}"} with
    \code{row.names} matching the \code{tile} names in \code{stgrid} and
    having the same \code{proj4string} as \code{events} and \code{W}. This is necessary
    to sample the spatial location of events generated by the endemic component.
  }
  \item{beta0,beta,gamma,siafpars,tiafpars}{
    these are the parameter subvectors of the \code{twinstim}.
    \code{beta} and \code{gamma} must be given in the 
    same order as they appear in \code{endemic} and \code{epidemic},
    respectively. \code{beta0} is either a single endemic intercept or a
    vector of type-specific endemic intercepts in the same order as in
    \code{qmatrix}.
  }
  \item{epilink}{
    a character string determining the link function to be used for the
    \code{epidemic} linear predictor of event marks. By default, the
    log-link is used. The experimental alternative is
    \code{epilink = "identity"}. Note that the identity link does not
    guarantee the force of infection to be positive. If this leads to a
    negative total intensity (endemic + epidemic), the point process is
    not well defined and simulation cannot proceed.
  }
  \item{t0}{
    \code{events} having occurred during (-Inf;\code{t0}] are regarded
    as part of the prehistory \eqn{H_0} of the process.
    For \code{simEpidataCS}, by default and also if \code{t0=NULL},
    the beginning of \code{stgrid} is used as \code{t0}.
    For the \code{simulate.twinstim} method, \code{NULL} means to use
    the fitted time range of the \code{"twinstim"} \code{object}.
  }
  \item{T, nEvents}{
    simulate a maximum of \code{nEvents} events up to time \code{T},
    then stop. For \code{simEpidataCS}, by default, and also if
    \code{T=NULL}, \code{T} equals the last stop time in \code{stgrid}
    (it cannot be greater) and \code{nEvents} is bounded above by 10000.
    For the \code{simulate.twinstim} method, \code{T=NULL} means to use
    the same same time range as for the fitting of the \code{"twinstim"}
    \code{object}. 
  } 
  \item{W}{
    see \code{\link{as.epidataCS}}. When simulating from
    \code{twinstim}-fits, \code{W} is by default taken from the original
    \code{data$W}. If specified as \code{NULL}, \code{W} is generated 
    automatically via \code{\link{unionSpatialPolygons}(tiles)}.
    However, since the result of such a polygon operation should always
    be verified, it is recommended to do that in advance.\cr
    It is important that \code{W} and \code{tiles} cover the same region:
    on the one hand direct offspring is sampled
    in the spatial influence region of the parent event, i.e., in the
    intersection of \code{W} and a circle of radius the \code{eps.s} of the
    parent event, after which the corresponding tile is determined by
    overlay with \code{tiles}. On the other hand endemic events are
    sampled from \code{tiles}.
  }
  \item{trace}{
    logical (or integer) indicating if (or how often) the current
    simulation status should be \code{cat}ed.  For the
    \code{simulate.twinstim} method, \code{trace} currently only applies
    to the first of the \code{nsim} simulations.
  }
  \item{.allocate}{
    number of rows (events) to initially allocate for the event history;
    defaults to 500.  Each time the simulated epidemic exceeds the
    allocated space, the event \code{data.frame} will be enlarged by
    \code{.allocate} rows.
  }
  \item{.skipChecks,.onlyEvents}{
    these logical arguments are not meant to be set by the user.
    They are used by the \code{simulate}-method for \code{"twinstim"} objects.
  }
  \item{object}{
    an object of class \code{"\link{twinstim}"}.
  }
  \item{nsim}{
    number of epidemics (i.e. spatio-temporal point patterns inheriting
    from class \code{"epidataCS"}) to simulate.  Defaults to 1 when the
    result is a simple object inheriting from class
    \code{"simEpidataCS"} (as if \code{simEpidataCS} would have been
    called directly).  If \code{nsim > 1}, the result 
    will be a list the structure of which depends on the argument
    \code{simplify}. 
  }
  \item{seed}{
    an object specifying how the random number generator should be
    initialized for simulation (via \code{\link{set.seed}}). The
    initial state will also be stored as an attribute \code{"seed"} of
    the result. The original state of the \code{\link{.Random.seed}}
    will be restored at the end of the simulation.
    By default (\code{NULL}), neither initialization nor recovery will
    be done.
    This behaviour is copied from the \code{\link{simulate}.lm} method.
  }
  \item{data}{
    an object of class \code{"epidataCS"}, usually the one to which the
    \code{"twinstim"} \code{object} was fitted. It carries
    the \code{stgrid} of the endemic component, but also
    \code{events} for use as the prehistory, and defaults for
    \code{rmarks} and \code{nCircle2Poly}.
  }
  \item{newcoef}{
    an optional named numeric vector of (a subset of) parameters to
    replace the original point estimates in \code{coef(object)}.
    Elements which do not match any model parameter by name are silently
    ignored. The \code{newcoef}s may also be supplied in a list
    following the same conventions as for the \code{start} argument in
    \code{\link{twinstim}}.
  }
  \item{simplify}{
    logical. It is strongly recommended to set \code{simplify = TRUE}
    (default) if \code{nsim} is large. This saves space and computation time,
    because for each simulated epidemic only the \code{events} component is
    saved. All other components, which do not vary between simulations,
    are only stored from the first run. In this case, the runtime of each
    simulation is stored as an attribute \code{"runtime"} to each simulated
    \code{events}. See also the \dQuote{Value} section below.
  }
  \item{control.siaf}{see \code{\link{twinstim}}.}
  \item{nCircle2Poly}{see \code{\link{as.epidataCS}}. For
    \code{simulate.twinstim}, \code{NULL} means to use the same value as
    for \code{data}.}
  \item{gmax}{
    maximum value the temporal interaction function
    \code{tiaf$g} can attain. If \code{NULL}, then it is assumed as the
    maximum value of the type-specific values at 0, i.e.
    \code{max(tiaf$g(rep.int(0,nTypes), tiafpars, 1:nTypes))}.
  }
  \item{\dots}{unused (arguments of the generic).}
}

\value{
  The function \code{simEpidataCS} returns a simulated epidemic of class
  \code{"simEpidataCS"}, which enhances the class
  \code{"epidataCS"} by the following additional components known from
  objects of class \code{"\link{twinstim}"}:
  \code{bbox}, \code{timeRange}, \code{formula}, \code{coefficients},
  \code{npars}, \code{control.siaf}, \code{call}, \code{runtime}.
  It has corresponding \code{\link{coeflist}},
  \code{\link[=residuals.simEpidataCS]{residuals}},
  \code{\link[=R0.simEpidataCS]{R0}}, and
  \code{\link[=intensityplot.simEpidataCS]{intensityplot}} methods.

  The \code{simulate.twinstim} method has some additional
  \emph{attributes} set on its result:
  \code{call}, \code{seed}, and \code{runtime}.
  If \code{nsim > 1}, it returns an object of class
  \code{"simEpidataCSlist"}, the form of which depends on the value of
  \code{simplify} (which is stored as an attribute \code{simplified}):
  if \code{simplify = FALSE}, then the return value is
  just a list of sequential simulations, each of class
  \code{"simEpidataCS"}. However, if \code{simplify = TRUE}, then the
  sequential simulations share all components but the simulated
  \code{events}, i.e. the result is a list with the same components as
  a single object of class \code{"simEpidataCS"}, but with \code{events}
  replaced by an \code{eventsList} containing the \code{events} returned
  by each of the simulations.

  The \code{stgrid} component of the returned \code{"simEpidataCS"}
  will be truncated to the actual end of the simulation, which might
  be \eqn{<T}, if the upper bound \code{nEvents} is reached during
  simulation.
  
  CAVE: Currently, \code{simplify=TRUE} in \code{simulate.twinstim}
  ignores that multiple simulated epidemics
  (\code{nsim > 1}) may have different \code{stgrid} 
  time ranges. In a \code{"simEpidataCSlist"}, the \code{stgrid} shared
  by all of the simulated epidemics is just the \code{stgrid}
  returned by the \emph{first} simulation.
}

\note{
  The more detailed the polygons in \code{tiles} are the slower is
  the algorithm. You are advised to sacrifice some shape
  details for speed by reducing the polygon complexity,
  for example via the \command{mapshaper} JavaScript library wrapped by
  the R package \CRANpkg{rmapshaper}, or via
  \code{\link[spatstat.geom]{simplify.owin}}.
}

\references{
  Douglas, D. H. and Peucker, T. K. (1973):
  Algorithms for the reduction of the number of points required to
  represent a digitized line or its caricature.
  \emph{Cartographica: The International Journal for Geographic
  Information and Geovisualization}, \bold{10}, 112-122

  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}
}

\author{ 
  Sebastian Meyer, with contributions by Michael H\enc{ö}{oe}hle
}

\seealso{
The function \code{\link{simEndemicEvents}} is a faster alternative
for endemic-only models, only returning a
\code{"\linkSPclass{SpatialPointsDataFrame}"} of simulated events.
  
The \code{\link{plot.epidataCS}} and \code{\link{animate.epidataCS}}
methods for plotting and animating continuous-space epidemic data,
respectively, also work for simulated epidemics (by inheritance),
and \code{\link{twinstim}} can be used to fit
spatio-temporal conditional intensity models also to simulated data.
}

\examples{
data("imdepi", "imdepifit")

## simulation needs the district polygons ('tiles')
load(system.file("shapes", "districtsD.RData", package="surveillance"))
if (surveillance.options("allExamples")) {
  plot(districtsD)
  plot(stateD, add=TRUE, border=2, lwd=2)
}

### simulate from the fitted model, only over the first 90 days (for speed)
nsim <- 10
\dontshow{if (!interactive()) nsim <- 2}% faster checks
sims <- simulate(imdepifit, nsim=nsim, seed=1, data=imdepi, T=90,
                 tiles=districtsD, simplify=TRUE)
sims

## plot event times for a random selection of 4 simulations
plot(sims, aggregate="time")

## extract the second realization -> object of class "simEpidataCS"
sim2 <- sims[[2]]
summary(sim2)
plot(sim2, aggregate="space")

## extract _cumulative_ number of events (simulated vs. observed)
getcsums <- function (events) {
    tapply(events$time, events@data["type"],
           function (t) cumsum(table(t)), simplify=FALSE)
}
csums_observed <- getcsums(imdepi$events)
csums_simulated <- lapply(sims$eventsList, getcsums)

## plot it
plotcsums <- function (csums, ...) {
    mapply(function (csum, ...) lines(as.numeric(names(csum)), csum, ...),
           csums, ...)
    invisible()
}
plot(c(0,90), c(0,35), type="n", xlab="Time [days]",
     ylab="Cumulative number of cases")
plotcsums(csums_observed, col=c(2,4), lwd=3)
legend("topleft", legend=levels(imdepi$events$type), col=c(2,4), lwd=1)
invisible(lapply(csums_simulated, plotcsums,
                 col=adjustcolor(c(2,4), alpha.f=0.5)))


\dontrun{
### Experimental code to generate 'nsim' simulations of 'nm2add' months
### beyond the observed time period:
nm2add <- 24
nsim <- 5
### The events still infective by the end of imdepi$stgrid will be used
### as the prehistory for the continued process.

origT <- tail(imdepi$stgrid$stop, 1)
## extend the 'stgrid' by replicating the last block 'nm2add' times
## (i.e., holding "popdensity" constant)
stgridext <- local({
    gLast <- subset(imdepi$stgrid, BLOCK == max(BLOCK))
    gAdd <- gLast[rep(1:nrow(gLast), nm2add),]; rownames(gAdd) <- NULL
    newstart <- seq(origT, by=30, length.out=nm2add)
    newstop <- c(newstart[-1], max(newstart) + 30)
    gAdd$start <- rep(newstart, each=nlevels(gAdd$tile))
    gAdd$stop <- rep(newstop, each=nlevels(gAdd$tile))
    rbind(imdepi$stgrid, gAdd, make.row.names = FALSE)[,-1]
})
## create an updated "epidataCS" with the time-extended 'stgrid'
imdepiext <- update(imdepi, stgrid = stgridext)
newT <- tail(imdepiext$stgrid$stop, 1)

## simulate beyond the original period
simsext <- simulate(imdepifit, nsim=nsim, seed=1, t0=origT, T=newT,
                    data=imdepiext, tiles=districtsD, simplify=TRUE)

## Aside to understand the note from checking events and tiles:
# marks(imdepi)["636",]  # tile 09662 is attributed to this event, but:
# plot(districtsD[c("09678","09662"),], border=1:2, lwd=2, axes=TRUE)
# points(imdepi$events["636",])
## this mismatch is due to polygon simplification

## plot the observed and simulated event numbers over time
plot(imdepiext, breaks=c(unique(imdepi$stgrid$start),origT),
     cumulative=list(maxat=330))
for (i in seq_along(simsext$eventsList))
    plot(simsext[[i]], add=TRUE, legend.types=FALSE,
         breaks=c(unique(simsext$stgrid$start),newT),
         subset=!is.na(source),  # have to exclude the events of the prehistory
         cumulative=list(offset=c(table(imdepi$events$type)), maxat=330, axis=FALSE),
         border=NA, density=0)  # no histogram
abline(v=origT, lty=2, lwd=2)

}
}

\keyword{datagen}
\keyword{models}
