\name{sts-class}
\docType{class}
\alias{sts}
\alias{sts-class}

% methods to access and replace slots
\alias{alarms,sts-method}
\alias{alarms<-,sts-method}
\alias{upperbound,sts-method}
\alias{upperbound<-,sts-method}
\alias{control,sts-method}
\alias{control<-,sts-method}
\alias{epoch,sts-method}
\alias{epoch<-,sts-method}
\alias{observed,sts-method}
\alias{observed<-,sts-method}
\alias{population,sts-method}
\alias{population<-,sts-method}
\alias{multinomialTS,sts-method}
\alias{multinomialTS<-,sts-method}
\alias{neighbourhood,sts-method}
\alias{neighbourhood<-,sts-method}
% other access methods
\alias{dim,sts-method}
\alias{dimnames,sts-method}
\alias{epochInYear}
\alias{epochInYear,sts-method}
\alias{year}
\alias{year,sts-method}
% constructor and conversion methods
\alias{initialize,sts-method}
\alias{as.data.frame,sts-method}
\alias{coerce,sts,ts-method}
\alias{coerce,ts,sts-method}
\alias{as.xts.sts}

\encoding{latin1}

\title{Class \code{"sts"} -- surveillance time series}

\description{

  This is a lightweight S4 class to implement (multivariate) time
  series of counts, typically from public health surveillance.
  For areal time series, the class can also capture the spatial layout
  of the regions, where the data originate from.
  The constructor function \code{sts} can be used to setup an
  \code{"sts"} object. Conversion of simple time-series objects
  (of class \code{"\link{ts}"}) is also possible.
  The slots of the \code{"sts"} class and available methods are
  described below.

}

\usage{
sts(observed, start = c(2000, 1), frequency = 52,
    population = NULL, ...)
}

\arguments{
  \item{observed}{a vector (for a single time series) or matrix (one
    time series per column) of counts. A purely numeric data frame will
    also do (transformed via \code{as.matrix}). This argument sets the
    \code{observed} slot, which is the core element of the resulting
    \code{"sts"} object. It determines the dimensions and colnames for
    several other slots. The columns (\dQuote{units}) typically
    correspond to different regions, diseases, or age groups.}
  \item{start,frequency}{basic characteristics of the time series data
    just like for simple \code{"\link{ts}"} objects. The (historical)
    default values correspond to weekly data starting in the first week
    of 2000.}
  \item{population}{a vector of length the number of columns in
    \code{observed} or a matrix of the same dimension as
    \code{observed}. Especially for multivariate time series, the
    population numbers (or fractions) underlying the counts in each unit
    are relevant for visualization and statistical inference. The
    \code{population} argument is an alias for the corresponding slot
    \code{populationFrac}. The default \code{NULL} value sets
    equal population fractions across all units.}
  \item{\dots}{further named arguments with names corresponding to slot
    names (see the list below). For instance, the \code{epoch} slot is
    used to store the observation time, either as an integer sequence
    (default) or as the numeric representation of \code{Date}s (if also
    \code{epochAsDate=TRUE}). In the public health surveillance context,
    the \code{state} slot is used to indicate outbreaks (default:
    \code{FALSE} for all observations). For areal time series data, the
    \code{map} and \code{neighbourhood} slots are used to store the
    spatial structure of the observation region.}
}

\section{Slots}{
  \describe{
    \item{\code{epoch}:}{Object of class \code{numeric} or specifying
      the time of observation. In old versions of the package this used to
      be the week numbers. However, depending on
      the \code{freq} argument, it can now be day or month as
      well. Furthermore, if \code{epochAsDate=TRUE} then 
      it is the \code{as.numeric} representation of \code{Date} objects
      giving the exact date of the observation. Note: This slot used to be
      called \code{week} in earlier versions of the package, but has now been
      renamed to reflect the greater 
      flexibility in the choice of observation time.}
    \item{\code{freq}:}{If weekly data \code{freq} corresponds to 52, in
      case of monthly data \code{freq} is 12.}
    \item{\code{start}:}{vector of length two denoting the year and the
      sample number (week, month, etc.) of the first observation}
    \item{\code{observed}:}{A matrix of size \code{length(epoch)} times the
      number of regions containing the weekly/monthly number of counts in
      each region. The colnames of the matrix should match the ID values of
      the shapes in the \code{map} slot.}
    \item{\code{state}:}{Matrix with the same dimension as \code{observed}
      containing booleans whether at the specific time point there was an
      outbreak in the region}
    \item{\code{alarm}:}{Matrix with the same dimension as
      \code{observed} specifying whether an outbreak detection algorithm
      declared a specific time point in the region as having an alarm.}
    \item{\code{upperbound}:}{Matrix with upper bound values }
    \item{\code{neighbourhood}:}{Symmetric matrix of size
      \eqn{(number of regions)^2} describing the neighbourhood structure. It
      may either be a binary adjacency matrix or contain neighbourhood orders
      (see the Examples for how to infer the latter from the \code{map}).}
    \item{\code{populationFrac}:}{A \code{matrix} of population 
      fractions or absolute numbers (see \code{multinomialTS} below)
      with dimensions \code{dim(observed)}.}
    \item{\code{map}:}{Object of class \code{SpatialPolygonsDataFrame}
      providing a shape of the areas which are monitored. }
    \item{\code{control}:}{Object of class \code{list}, this is a
      rather free data type to be returned by the surveillance algorithms. }
    \item{\code{epochAsDate}:}{Object of class \code{"logical"} stating
      whether to use a ISO 8601 representation of the \code{epoch} slot using
      the \code{Date} class (\code{epochAsDate=TRUE}) or just to
      interpret the epochs as numerics (\code{epochAsDate=FALSE}).}
    \item{\code{multinomialTS}:}{Object of class \code{"logical"}
      stating whether to interpret the object as \code{observed} out of
      \code{population}, i.e. a multinomial interpretation instead of a
      count interpretation.}
  }
}

\section{Methods}{
  \describe{ 
    \item{dim}{\code{signature(x = "sts")}: extract matrix dimensions of
      \code{observed}.}
    \item{dimnames}{\code{signature(x="sts")}: extract the \code{\link{dimnames}} of the \code{observed} matrix.}
    \item{observed}{\code{signature(x = "sts")}: extract the \code{observed}
      slot of an \code{sts} object.}
    \item{population}{\code{signature(x = "sts")}: extract the \code{populationFrac}
      slot of an \code{sts} object.}
    \item{multinomialTS}{\code{signature(x = "sts")}: extract the
      \code{multinomialTS} slot of an \code{sts} object.}
    \item{neighbourhood}{\code{signature(x = "sts")}: extract the
	  \code{neighbourhood} slot of an \code{sts} object.}
    \item{alarms}{\code{signature(x = "sts")}: extract the \code{alarm}
      slot of an \code{sts} object.}
    \item{upperbound}{\code{signature(x = "sts")}: extract the \code{upperbound}
      slot of an \code{sts} object.}
    \item{control}{\code{signature(x = "sts")}: extract the \code{control}
      slot of an \code{sts} object.}
    \item{epoch}{\code{signature(x = "sts")}: extract the \code{epoch}
      slot of an \code{sts} object. If ISO dates are used then the
      returned object is of class \code{Date}.}
    \item{epochInYear}{\code{signature(x = "sts")}: Returns the epoch
      number within the year of the \code{epoch} slot.}
    \item{initialize}{\code{signature(x="sts")}: the internal function
      \code{init.sts} is called, which assigns all slots.  }
    \item{[}{subsetting \code{"sts"} objects, see \code{\link{[,sts-method}}.}
    \item{aggregate}{\code{signature(x="sts")}: see
      \code{\link{aggregate,sts-method}}}
    \item{year}{\code{signature(x = "sts")}:
      extracts the corresponding year of each observation of \code{x}}
    \item{as.data.frame}{\code{signature(x = "sts")}:
      converts the \code{observed}, \code{epoch}, \code{state} and
      \code{alarm} slots of \code{x} into a data frame with column names
      matching the colnames of the respective slots. Useful when one wants
      to fit a model based on the object}
    \item{plot}{\code{signature(x="sts",y="missing")}: this method is
      the entry point to a collection of plot variants. It is also
      the successor of the \code{\link{plot.disProg}} and 
      \code{\link{plot.survRes}} functions. The type of plot is
      specified using a formula \code{type}.
      See \code{\link{stsplot}} for details.
    }
    \item{animate}{see \code{\link{animate.sts}}.}
    \item{coerce}{\code{signature(from="sts", to="ts")} (can also be
      called via the S3 generic \code{as.ts}) and the reverse
      \code{signature(from="ts", to="sts")}.}
    \item{toLatex}{see \code{\link{toLatex,sts-method}}.}
  }
}

\author{Michael \enc{Höhle}{Hoehle} and Sebastian Meyer}

\examples{
## A typical dataset with weekly counts of measles from several districts
data("measlesWeserEms")
measlesWeserEms

## reconstruct data("measlesWeserEms") from its components
counts <- observed(measlesWeserEms)
map <- measlesWeserEms@map
populationFrac <- population(measlesWeserEms)
weserems_nbOrder <- neighbourhood(measlesWeserEms)
## orders of adjacency can also be determined from the map
if (requireNamespace("spdep")) {
    stopifnot(identical(weserems_nbOrder,
                        nbOrder(poly2adjmat(map), maxlag = 10)))
}
mymeasles <- sts(counts, start = c(2001, 1), frequency = 52,
                 population = populationFrac,
                 neighbourhood = weserems_nbOrder, map = map)
stopifnot(identical(mymeasles, measlesWeserEms))

## convert ts/mts object to sts
z <- ts(matrix(rpois(300,10), 100, 3), start = c(1961, 1), frequency = 12)
z.sts <- as(z, "sts")
plot(z.sts)

## conversion of "sts" objects to the quasi-standard "xts" class
if (requireNamespace("xts")) {
    z.xts <- as.xts.sts(z.sts)
    plot(z.xts)
}
}

\keyword{classes}
