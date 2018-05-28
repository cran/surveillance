\encoding{latin1}
\name{epidataCS}
\alias{epidataCS}
\alias{as.epidataCS}
\alias{print.epidataCS}
\alias{nobs.epidataCS}
\alias{head.epidataCS}
\alias{tail.epidataCS}
\alias{[.epidataCS}
\alias{subset.epidataCS}
\alias{marks.epidataCS}
\alias{summary.epidataCS}
\alias{print.summary.epidataCS}
\alias{as.stepfun.epidataCS}
\alias{getSourceDists}
\alias{coerce,epidataCS,SpatialPointsDataFrame-method}

\title{
  Continuous Space-Time Marked Point Patterns with Grid-Based Covariates
}

\description{
  Data structure for \strong{c}ontinuous \strong{s}patio-temporal event
  data, e.g. individual case reports of an infectious disease.
  Apart from the actual \code{events}, the class simultaneously
  holds a spatio-temporal grid of endemic covariates (similar to
  disease mapping) and a representation of the observation region.

  The \code{"epidataCS"} class is the basis for fitting 
  spatio-temporal endemic-epidemic intensity models with the function
  \code{\link{twinstim}} (Meyer et al., 2012).
  The implementation is described in Meyer et al. (2017, Section 3),
  see \code{vignette("twinstim")}.
}

\usage{
as.epidataCS(events, stgrid, W, qmatrix = diag(nTypes),
             nCircle2Poly = 32L, T = NULL,
             clipper = c("polyclip", "rgeos"), verbose = interactive())

\method{print}{epidataCS}(x, n = 6L, digits = getOption("digits"), ...)

\method{nobs}{epidataCS}(object, ...)
\method{head}{epidataCS}(x, n = 6L, ...)
\method{tail}{epidataCS}(x, n = 6L, ...)
\method{[}{epidataCS}(x, i, j, ..., drop = TRUE)
\method{subset}{epidataCS}(x, subset, select, drop = TRUE, ...)

\method{marks}{epidataCS}(x, coords = TRUE, ...)

\method{summary}{epidataCS}(object, ...)
\method{print}{summary.epidataCS}(x, ...)

\method{as.stepfun}{epidataCS}(x, ...)

getSourceDists(object, dimension = c("space", "time"))
}

\arguments{
  \item{events}{
    a \code{"\linkS4class{SpatialPointsDataFrame}"} of cases with the
    following obligatory columns (in the \code{events@data}
    \code{data.frame}): 
    \describe{
      \item{time}{time point of event. Will be converted to a numeric
        variable by \code{as.numeric}. There should be no concurrent
        events (but see \code{\link{untie}} for an ex post adjustment)
        and there cannot be events beyond \code{stgrid}
        (i.e., \code{time<=T} is required). Events at or before time
        \eqn{t_0} = \code{min(stgrid$start)} are allowed and form the
        prehistory of the process.}
      \item{tile}{the spatial region (tile) where the event is located.
      This links to the tiles of \code{stgrid}.}
      \item{type}{optional type of event in a marked \code{twinstim}
	model. Will be converted to a factor variable dropping unused
        levels. If missing, all events will be attribute the single type
        \code{"1"}.}
      \item{eps.t}{maximum \emph{temporal} influence radius (e.g. length of
	infectious period, time to culling, etc.); must be positive and
        may be \code{Inf}.} 
      \item{eps.s}{maximum \emph{spatial} influence radius (e.g. 100 [km]);
        must be positive and may be \code{Inf}. A compact influence
        region mainly has computational advantages, but might also be
        plausible for specific applications.}
    }
    The \code{data.frame} may contain columns with further marks of
    the events, e.g. sex, age of infected individuals, which may
    be used as epidemic covariates influencing infectiousness.
    Note that some auxiliary columns will be added at conversion
    whose names are reserved: \code{".obsInfLength"},
    \code{".bdist"}, \code{".influenceRegion"}, and \code{".sources"},
    as well as \code{"start"}, \code{"BLOCK"}, and all endemic
    covariates' names from \code{stgrid}.
  }

  \item{stgrid}{
    a \code{\link{data.frame}} describing endemic covariates on a full
    spatio-temporal region x interval grid (e.g., district x week),
    which is a decomposition of the observation region \code{W} and
    period \eqn{t_0,T}. This means that for every combination of spatial
    region and time interval there must be exactly one row in this
    \code{data.frame}, that the union of the spatial tiles equals
    \code{W}, the union of the time intervals equals \eqn{t_0,T}, and
    that regions (and intervals) are non-overlapping.
    There are the following obligatory columns: 
    \describe{
      \item{tile}{ID of the spatial region (e.g., district ID). It will
        be converted to a factor variable (dropping unused levels if it
        already was one).}
      \item{start, stop}{columns describing the consecutive temporal
        intervals (converted to numeric variables by \code{as.numeric}).
        The \code{start} time of an interval must be equal to the
        \code{stop} time of the previous interval. The \code{stop} column may
	be missing, in which case it will be auto-generated from the set
	of \code{start} values and \code{T}.}
      \item{area}{area of the spatial region (\code{tile}).
	Be aware that the unit of this area (e.g., square km) must be consistent
	with the units of \code{W} and \code{events} (as specified in
	their \code{\link{proj4string}}s, if they have projected coordinates).}
    }
    The remaining columns are endemic covariates.
    Note that the column name \code{"BLOCK"} is reserved
    (a column which will be added automatically for indexing the time
    intervals of \code{stgrid}).
  }

  \item{W}{
    an object of class \code{"\linkS4class{SpatialPolygons}"}
    representing the observation region.
    It must have the same \code{proj4string} as \code{events}
    and all events must be within \code{W}.
    Prior simplification of \code{W} may considerably reduce the
    computational burden of likelihood evaluations in
    \code{\link{twinstim}} models with non-trivial spatial
    interaction functions (see the \dQuote{Note} section below).
  }

  \item{qmatrix}{
    a square indicator matrix (0/1 or \code{FALSE}/\code{TRUE}) for possible
    transmission between the event types. The matrix will be internally
    converted to \code{logical}. Defaults to an independent spread of the event
    types, i.e. the identity matrix. 
  }

  \item{nCircle2Poly}{
    accuracy (number of edges) of the polygonal approximation of a circle,
    see \code{\link{discpoly}}.
  }

  \item{T}{
    end of observation period (i.e. last \code{stop} time of
    \code{stgrid}). Must be specified if the start but not the stop
    times are supplied in \code{stgrid} (=> auto-generation of
    \code{stop} times).
  }

  \item{clipper}{polygon clipping engine to use for calculating the
    \code{.influenceRegion}s of events (see the Value section below).
    Default is the \CRANpkg{polyclip} package (called via
    \code{\link{intersect.owin}} from package \CRANpkg{spatstat}).
    In \pkg{surveillance} <= 1.6-0, package \pkg{gpclib} was used, which
    has a restrictive license. This is no longer supported.}

  \item{verbose}{logical indicating if status messages should be printed
  during input checking and \code{"epidataCS"} generation. The default
  is to do so in interactive \R sessions.}

  \item{x}{an object of class \code{"epidataCS"} or
    \code{"summary.epidataCS"}, respectively.}

  \item{n}{a single integer. If positive, the first (\code{head}, \code{print})
    / last (\code{tail}) \code{n} events are extracted. If negative,
    all but the \code{n} first/last events are extracted.}

  \item{digits}{minimum number of significant digits to be printed in
    values.}

  \item{i,j,drop}{
    arguments passed to the
    \code{\link[=[,SpatialPointsDataFrame-method]{[-method}} for
    \code{SpatialPointDataFrame}s for subsetting the \code{events} while
    retaining \code{stgrid} and \code{W}.\cr
    If \code{drop=TRUE} (the default), event types that completely
    disappear due to \code{i}-subsetting will be dropped, which reduces
    \code{qmatrix} and the factor levels of the \code{type} column.\cr
    By the \code{j} index, epidemic covariates can be removed from
    \code{events}.}

  \item{\dots}{unused (arguments of the generics) with a few exceptions:
    The \code{print} method for \code{"epidataCS"} passes
    \code{\dots} to the \code{\link{print.data.frame}} method, and the
    \code{print} method for \code{"summary.epidataCS"} passes additional
    arguments to \code{\link{print.table}}.}

  \item{subset, select}{arguments used to subset the \code{events} from
    an \code{"epidataCS"} object like in \code{\link{subset.data.frame}}.}

  \item{coords}{logical indicating if the data frame of event marks
    returned by \code{marks(x)} should have the event
    coordinates appended as last columns. This defaults to \code{TRUE}.}

  \item{object}{an object of class \code{"epidataCS"}.}

  \item{dimension}{the distances of all events to their potential source
    events can be computed in either the \code{"space"} or \code{"time"}
    dimension.}
}

\details{
  The function \code{as.epidataCS} is used to generate objects of class
  \code{"epidataCS"}, which is the data structure required for
  \code{\link{twinstim}} models.

  The \code{[}-method for class \code{"epidataCS"}
  ensures that the subsetted object will be valid, for instance, it
  updates the auxiliary list of potential transmission paths stored
  in the object. The \code{[}-method is used in
  \code{subset.epidataCS}, which is implemented similar to
  \code{\link{subset.data.frame}}.
  
  The \code{print} method for \code{"epidataCS"} prints some metadata
  of the epidemic, e.g., the observation period, the dimensions of the
  spatio-temporal grid, the types of events, and the total number of
  events. By default, it also prints the first \code{n = 6} rows of the
  \code{events}.
}

\value{
  An object of class \code{"epidataCS"} is a list containing the
  following components: 
  \item{events}{a \code{"\linkS4class{SpatialPointsDataFrame}"} (see the
    description of the argument). 
    The input \code{events} are checked for requirements and sorted
    chronologically. The columns are in the following
    order: obligatory event columns, event marks, the columns \code{BLOCK},
    \code{start} and endemic covariates copied from \code{stgrid},
    and finally, hidden auxiliary columns. 
    The added auxiliary columns are:
    \describe{
      \item{\code{.obsInfLength}}{observed length of the infectious period
	(possibly truncated at \code{T}), i.e., \code{pmin(T-time, eps.t)}.}
      \item{\code{.sources}}{a list of numeric vectors of potential sources of
        infection (wrt the interaction ranges eps.s and eps.t) for each
        event. Row numbers are used as index.}
      \item{\code{.bdist}}{minimal distance of the event locations to the
        polygonal boundary \code{W}.}
      \item{\code{.influenceRegion}}{a list of influence regions represented by
        objects of the \pkg{spatstat} class \code{"owin"}. For each
        event, this is the intersection of \code{W} with a (polygonal)
        circle of radius \code{eps.s} centered at the event's location,
        shifted such that the event location becomes the origin.
        The list has \code{nCircle2Poly} set as an attribute.}
      }
  }
  \item{stgrid}{a \code{data.frame} (see description of the argument).
    The spatio-temporal grid of endemic covariates is sorted by time
    interval (indexed by the added variable \code{BLOCK}) and region
    (\code{tile}). It is a full \code{BLOCK} x \code{tile} grid.}
  \item{W}{a \code{"\linkS4class{SpatialPolygons}"} object representing
    the observation region.}
  \item{qmatrix}{see the above description of the argument. The
    \code{\link{storage.mode}} of the indicator matrix is set to logical
    and the \code{dimnames} are set to the levels of the event types.}

  The \code{nobs}-method returns the number of events.

  The \code{head} and \code{tail} methods subset the epidemic data using
  the extraction method (\code{[}), i.e. they return an object of class
  \code{"epidataCS"}, which only contains (all but) the first/last
  \code{n} events.

  For the \code{"epidataCS"} class, the method of the generic function
  \code{\link[spatstat]{marks}} defined by the \pkg{spatstat} package
  returns a \code{data.frame} of the event marks (actually also
  including time and location of the events), disregarding endemic
  covariates and the auxiliary columns from the \code{events} component
  of the \code{"epidataCS"} object.

  The \code{summary} method (which has again a \code{print} method)
  returns a list of metadata, event data, the tables of tiles and types,
  a step function of the number of infectious individuals over time
  (\code{$counter}), i.e., the result of the
  \code{\link{as.stepfun}}-method for \code{"epidataCS"}, and the number
  of potential sources of transmission for each event (\code{$nSources})
  which is based on the given maximum interaction ranges \code{eps.t}
  and \code{eps.s}.
}

\note{
  Since the observation region \code{W} defines the integration domain
  in the point process likelihood,
  the more detailed the polygons of \code{W} are the longer it will
  take to fit a \code{\link{twinstim}}. You are advised to
  sacrifice some shape details for speed by reducing the polygon
  complexity, for example via \code{\link[rmapshaper]{ms_simplify}} from
  the \CRANpkg{rmapshaper} package.
  Alternative tools are provided by the packages \CRANpkg{maptools}
  (\code{\link[maptools]{thinnedSpatialPoly}}) and \CRANpkg{spatstat}
  (\code{\link[spatstat]{simplify.owin}}).
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}

  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
}

\author{
  Sebastian Meyer

  Contributions to this documentation by Michael
  H\enc{ö}{oe}hle and Mayeul Kauffmann.
}

\seealso{
  \code{vignette("twinstim")}.

  \code{\link{plot.epidataCS}} for plotting, and
  \code{\link{animate.epidataCS}} for the animation of such an epidemic.
  There is also an \code{\link[=update.epidataCS]{update}} method for the
  \code{"epidataCS"} class.

  To re-extract the \code{events} point pattern from \code{"epidataCS"},
  use \code{as(object, "SpatialPointsDataFrame")}.

  It is possible to convert an \code{"epidataCS"} point pattern to
  an \code{"\link{epidata}"} object (\code{\link{as.epidata.epidataCS}}),
  or to aggregate the events into an \code{"\linkS4class{sts}"} object
  (\code{\link{epidataCS2sts}}).
}

\examples{
## load "imdepi" example data (which is an object of class "epidataCS")
data("imdepi")

## print and summary
print(imdepi, n=5, digits=2)
print(s <- summary(imdepi))
plot(s$counter,  # same as 'as.stepfun(imdepi)'
     xlab = "Time [days]", ylab="Number of infectious individuals",
     main=paste("Time course of the number of infectious individuals",
                "assuming an infectious period of 30 days", sep="\n"))
plot(table(s$nSources), xlab="Number of \"close\" infective individuals",
     ylab="Number of events",
     main=paste("Distribution of the number of potential sources",
                "assuming an interaction range of 200 km and 30 days",
                sep="\n"))
## the summary object contains further information
str(s)

## a histogram of the spatial distances to potential source events
## (i.e., to events of the previous eps.t=30 days within eps.s=200 km)
sourceDists_space <- getSourceDists(imdepi, "space")
hist(sourceDists_space); rug(sourceDists_space)

## internal structure of an "epidataCS"-object
str(imdepi, max.level=4)
## see help("imdepi") for more info on the data set

## extraction methods subset the 'events' component
## (thereby taking care of the validity of the epidataCS object,
## for instance the hidden auxiliary column .sources)
imdepi[101:200,]
tail(imdepi, n=4)           # reduce the epidemic to the last 4 events
subset(imdepi, type=="B")   # only consider event type B

## see help("plot.epidataCS") for convenient plot-methods for "epidataCS"


###
### reconstruct the "imdepi" object
###

## observation region
load(system.file("shapes", "districtsD.RData", package="surveillance"),
     verbose = TRUE)
summary(stateD)

## extract point pattern of events from the "imdepi" data
data(imdepi)
events <- marks(imdepi)  # data frame with coordinate columns
coordinates(events) <- c("x", "y")  # promote to a "SpatialPointsDataFrame"
#proj4string(events) <- proj4string(stateD)
events@proj4string <- stateD@proj4string  # exact copy (avoid CRS reformatting)

## or, much simpler, use the corresponding coerce-method
\dontshow{
events@coords.nrs <- numeric(0L)
stopifnot(all.equal(as(imdepi, "SpatialPointsDataFrame"), events))
}
events <- as(imdepi, "SpatialPointsDataFrame")
summary(events)

## plot observation region with events
plot(stateD, axes=TRUE); title(xlab="x [km]", ylab="y [km]")
points(events, pch=unclass(events$type), cex=0.5, col=unclass(events$type))
legend("topright", legend=levels(events$type), title="Type", pch=1:2, col=1:2)

## space-time grid with endemic covariates
head(stgrid <- imdepi$stgrid[,-1])

## reconstruct the "imdepi" object from its components
myimdepi <- as.epidataCS(events = events, stgrid = stgrid,
                         W = stateD, qmatrix = diag(2), nCircle2Poly = 16)
\dontrun{
## This reconstructed object is equal to 'imdepi' as long as the internal
## structures of the embedded classes ("owin", "SpatialPolygons", ...), and
## the calculation of the influence regions by "polyclip" have not changed:
stopifnot(all.equal(imdepi, myimdepi))
}
}

\keyword{spatial}
\keyword{classes}
\keyword{manip}
