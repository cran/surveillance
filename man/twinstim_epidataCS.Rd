\encoding{latin1}
\name{twinstim_epidataCS}
\alias{epidataCS}
\alias{as.epidataCS}
\alias{print.epidataCS}
\alias{nobs.epidataCS}
\alias{head.epidataCS}
\alias{tail.epidataCS}
\alias{[.epidataCS}
\alias{subset.epidataCS}
\alias{marks}	       % the generic 'marks()' is imported from spatstat and
                       % then exported (unchanged) -> documentation required
\alias{marks.ppp}               % dito
\alias{markformat.default}      % dito
\alias{marks.epidataCS}
\alias{summary.epidataCS}
\alias{print.summary.epidataCS}

\title{
  Class for Representing Continuous Space-Time Point Process Data
}

\description{
  Data structure for (c)ontinuous (s)patio-temporal event data, typically
  cases of an infectious disease. The data structure simultaneously
  contains a spatio-temporal grid of endemic covariates and a
  representation of the observation region.

  The \code{"epidataCS"} class is the basis for fitting 
  spatio-temporal epidemic intensity models with the function
  \code{\link{twinstim}}.
}

\usage{
as.epidataCS(events, stgrid, W, qmatrix = diag(nTypes),
             nCircle2Poly = 32L, T = NULL)

\method{print}{epidataCS}(x, n = 6L, digits = getOption("digits"), ...)

\method{nobs}{epidataCS}(object, ...)
\method{head}{epidataCS}(x, n = 6L, ...)
\method{tail}{epidataCS}(x, n = 6L, ...)
\method{[}{epidataCS}(x, i, j, drop = FALSE)
\method{subset}{epidataCS}(x, subset, select, drop = FALSE, ...)

\method{marks}{epidataCS}(x, coords = TRUE, ...)

\method{summary}{epidataCS}(object, ...)
\method{print}{summary.epidataCS}(x, ...)
}

\arguments{
  \item{events}{
    \code{"\linkS4class{SpatialPointsDataFrame}"} of cases with the
    following obligatory columns (in the \code{events@data}
    \code{data.frame}): 
    \describe{
      \item{time}{time point of event. Will be converted to a numeric
        variable by \code{as.numeric}. There should be no concurrent
        events and the event times must be covered by \code{stgrid},
        i.e. belong to the time interval \eqn{(t_0,T]}, where \eqn{t_0}
        is \code{min(stgrid$start)} and \code{T} is described below.}
      \item{tile}{reference to the spatial region (tile) in \code{stgrid},
	where the event is located.}
      \item{type}{optional type of event in a marked \code{twinstim}
	model. Will be converted to a factor variable dropping unused
        levels. If missing, all events will be attribute the single type
        \code{"1"}.}
      \item{eps.t}{maximum \emph{temporal} influence radius (e.g. length of
	infectious period, time to culling, etc.); must be positive and
        may be \code{Inf}.} 
      \item{eps.s}{maximum \emph{spatial} influence radius (e.g. 100 [km]);
        must be positive and may be \code{Inf}.}
    }
    The \code{data.frame} may contain columns with further marks of
    the events, e.g. sex, age of infected individuals, which may
    be used as epidemic covariates influencing infectiousness.
    Note that some auxiliary columns will be added at conversion
    whose names are reserved: \code{"ID"}, \code{".obsInfLength"},
    \code{".bdist"}, \code{".influenceRegion"}, and \code{".sources"},
    as well as \code{"start"}, \code{"BLOCK"}, and all endemic
    covariates' names from \code{stgrid}.
  }

  \item{stgrid}{
    \code{data.frame} describing endemic covariates on a full
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
    an object of class \code{"\linkS4class{SpatialPolygons}"} indicating the
    observation region. Must have the same \code{proj4string} as \code{events}
    and all events must be within \code{W}.
    The function \code{\link{simplify.owin}} 
    from the package \pkg{spatstat} may be useful if polygonal operations take
    too long or memory is limited (see also the \dQuote{Note} section below).
  }

  \item{qmatrix}{
    square indicator matrix (0/1 or \code{FALSE}/\code{TRUE}) for possible
    transmission between the event types. The matrix will be internally
    converted to \code{logical}. Defaults to an independent spread of the event
    types, i.e. the identity matrix. 
  }

  \item{nCircle2Poly}{
    accuracy (number of edges) of the polygonal approximation of a circle.
  }

  \item{T}{
    end of observation period (i.e. last \code{stop} time of
    \code{stgrid}). Must be specified if the start but not the stop
    times are supplied in \code{stgrid} (=> auto-generation of
    \code{stop} times).
  }

  \item{x}{an object of class \code{"epidataCS"} or
    \code{"summary.epidataCS"}, respectively.}

  \item{n}{a single integer. If positive, the first (\code{head}, \code{print})
    / last (\code{tail}) \code{n} events are extracted. If negative,
    all but the \code{n} first/last events are extracted.
  }

  \item{digits}{minimum number of significant digits to be printed in
    values.}

  \item{i,j}{
    arguments passed to the
    \code{\link[=[,SpatialPointsDataFrame-method]{[-method}} for
    \code{SpatialPointDataFrame}s for subsetting the \code{events} (not
    \code{stgrid} or \code{W}).
    Only epidemic covariates can be removed from \code{events} by the
    \code{j} index. The other components of \code{x} (\code{stgrid},
    \ldots) are retained.}

  \item{drop, \dots}{unused (arguments of the generics).
    However, the \code{print} method for class \code{"epidataCS"} passes
    \dots to the \code{\link{print.data.frame}} method.}

  \item{subset, select}{arguments used to subset the \code{events} from
    an \code{"epidataCS"} object like in \code{\link{subset.data.frame}}.}

  \item{coords}{logical indicating if the data frame of event marks
    returned by \code{marks.epidataCS} should also include the event
    coordinates. This defaults to \code{TRUE}.}

  \item{object}{an object of class \code{"epidataCS"}.}
}

\details{
  The function \code{as.epidataCS} is used to generate objects of class
  \code{"epidataCS"}, which is the data structure that
  \code{\link{twinstim}} models work on.

  The extraction method for class \code{"epidataCS"}
  ensures that the subsetted object will be valid, for instance, it
  updates the auxiliary list of potential transmission paths stored
  in the object. This \code{[}-method is also the basis for the
  \code{subset.epidataCS}-method, which is implemented similar to the
  \code{\link{subset.data.frame}}-method.
  
  The \code{print} method for \code{"epidataCS"} prints some metadata
  of the epidemic, e.g., the observation period, the dimensions of the
  spatio-temporal grid, the types of events, and the total number of
  events. By default, it also prints the first \code{n = 6} rows of the
  \code{events}.
}

\value{
  An object of class \code{"epidataCS"} is a list containing the
  following components: 
  \item{events}{see the description of the argument.
    The input \code{events} are checked for requirements, sorted
    chronologically and attributed an ID. The columns are in the following
    order: ID, obligatory event columns, event marks, the columns \code{BLOCK},
    \code{start} and endemic covariates copied from \code{stgrid},
    and finally, hidden auxiliary columns. 
    The added auxiliary columns are:
    \describe{
      \item{.obsInfLength}{observed length of the infectious period
	(being part [0,T]), i.e. \code{pmin(T-time, eps.t)}.}
      \item{.sources}{a list of numeric vectors of potential sources of
        infection (wrt the interaction ranges eps.s and eps.t) for each
        event. Currently row numbers are used as index, not the ID's.}
      \item{.bdist}{minimal distance of the event locations to the
        polygonal boundary \code{W}.}
      \item{.influenceRegion}{a list of influence regions represented by
        objects of the \pkg{spatstat} class \code{"owin"}. For each
        event, this is the intersection of \code{W} with a (polygonal)
        circle of radius \code{eps.s} centered at the event location.
        The list has \code{nCircle2Poly} set as an attribute.}
      }
  }
  \item{stgrid}{see the description of the argument. The spatio-temporal
    grid of endemic covariates is sorted by time interval (indexed by the
    added variable \code{BLOCK}) and region (\code{tile}). It is a
    full \code{BLOCK} x \code{tile} grid.}
  \item{W}{observation region as \code{SpatialPolygons} (see the
    argument description).}
  \item{qmatrix}{see the description of the argument. The
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
  (\code{counter}), and the number of
  potential sources of transmission for each event (\code{nSources})
  which is based on the given maximum interaction ranges \code{eps.t}
  and \code{eps.s}.
}

\note{
  The more detailed the observation region \code{W} is the slower the
  operations will be. Often it can be advantageous to sacrifice some detail
  for speed by reducing polygon complexity using, e.g., the
  Douglas and Peucker (1973) reduction method available at
  \url{http://MapShaper.org} (Harrower and Bloch, 2006) or as function
  \code{dp()} in the (slightly outdated) package \pkg{shapefiles}, or by
  passing by \pkg{spatstat}'s \code{\link[spatstat]{simplify.owin}} procedure.

  Note also that generation of \code{"epidataCS"} depends on the
  computation of polygon intersections via the \pkg{gpclib} package,
  which has a restricted license. In order to use \code{as.epidataCS},
  this license must be explicitly accepted by setting
  \code{\link{surveillance.options}(gpclib=TRUE)}.
}

\references{
  Douglas, D. H. and Peucker, T. K. (1973):
  Algorithms for the reduction of the number of points required to
  represent a digitized line or its caricature.
  \emph{Cartographica: The International Journal for Geographic
  Information and Geovisualization}, \bold{10}, 112-122.

  Harrower, M. and Bloch, M. (2006):
  MapShaper.org: A Map Generalization Web Service.
  \emph{IEEE Computer Graphics and Applications}, \bold{26}(4), 22-27.\cr
  DOI-Link: \url{http://dx.doi.org/10.1109/MCG.2006.85}

  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}

  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{ä}{ae}t
  M\enc{ü}{ue}nchen.\cr
  Available as \url{http://epub.ub.uni-muenchen.de/11703/}
}

\author{
  Sebastian Meyer with documentation contributions by Michael
  H\enc{ö}{oe}hle and Mayeul Kauffmann.
}

\seealso{
  \code{\link{plot.epidataCS}} for plotting, and
  \code{\link{animate.epidataCS}} for the animation of such an epidemic.
  There is also an \code{\link[=update.epidataCS]{update}} method for the
  \code{"epidataCS"} class.
  Models for \code{"epidataCS"} can be fitted with \code{\link{twinstim}}.
  It is also possible to convert the data to \code{\link{epidata}}
  objects (discrete space) for analysis with \code{\link{twinSIR}}
  (see \code{\link{as.epidata.epidataCS}}).
}

\examples{
## load "imdepi" example data (which is an object of class "epidataCS")
data(imdepi)

## take a look at the data object
print(imdepi, n=5, digits=2)
s <- summary(imdepi)
s
plot(s$counter, xlab = "Time [days]",
     ylab="Number of infectious individuals",
     main=paste("Time series of the number of infectious individuals",
                "assuming an infectious period of 30 days", sep="\n"))
plot(table(s$nSources), xlab="Number of \"close\" infective individuals",
     ylab="Number of events",
     main=paste("Distribution of the number of potential sources",
                "assuming an interaction range of 200 km and 30 days",
                sep="\n"))
## the summary object contains further information
str(s)

## see the help page on the 'imdepi' dataset for more insight

## extraction methods subset the 'events' component
## (thereby taking care of the validity of the epidataCS object,
## for instance the hidden auxiliary column .sources and qmatrix)
imdepi[101:200, -match("sex", names(imdepi$events))]
tail(imdepi, n=4)           # reduce the epidemic to the last 4 events
subset(imdepi, type=="B")   # only consider event type B


### now reconstruct the object from its components

## (Simplified) observation region (as SpatialPolygons)
W <- imdepi$W
summary(W)

## events
events <- marks(imdepi)[,-1]
coordinates(events) <- c("x", "y")    # promote to a SpatialPointsDataFrame
proj4string(events) <- proj4string(W) # ETRS89 projection (+units=km)
summary(events)

## plot W with events
plot(W, axes=TRUE)
title(xlab="x [km]", ylab="y [km]")
points(events, pch=unclass(events$type), cex=0.5, col=unclass(events$type))
legend("topright", legend=levels(events$type), title="Type",
       pch=1:2, col=1:2)
\dontrun{
    # similar to using the convenient plot-method for "epidataCS"
    plot(imdepi, "space")
}

## endemic covariates
stgrid <- imdepi$stgrid[,-1]

## indicator matrix for possible transmission between event types
qmatrix <- imdepi$qmatrix
qmatrix

## now construct an object of class "epidataCS" from these components
if (require("gpclib")) {
    oopt <- surveillance.options(gpclib=TRUE)
    myimdepi <- as.epidataCS(events, stgrid, W, qmatrix, nCircle2Poly=16)
    surveillance.options(oopt)
}

\dontshow{
## should be equal as long as the internal structures of "epidataCS" and
## of the embedded class "owin" do not change
if (!identical(imdepi, myimdepi))
    warning("@Sebastian: update data(imdepi)")
}

## take a look into the internal structure of an epidataCS-object
str(imdepi, max.level=4)
}

\keyword{spatial}
\keyword{classes}
\keyword{manip}
