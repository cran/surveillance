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
\alias{frequency,sts-method}
\alias{start,sts-method}
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
\alias{year}
\alias{year,sts-method}
\alias{epochInYear}
\alias{epochInYear,sts-method}
% conversion methods
\alias{as.data.frame.sts}
\alias{as.data.frame,sts-method}
\alias{as.ts.sts}
\alias{coerce,sts,ts-method}
\alias{coerce,ts,sts-method}
\alias{as.xts.sts}

\encoding{latin1}

\title{Class \code{"sts"} -- surveillance time series}

\description{

  This is a lightweight S4 class to implement (multivariate) time
  series of counts, typically from public health surveillance.
  The \code{"sts"} class supersedes the informal \code{"disProg"} class
  used in early versions of package \pkg{surveillance}. Converters are
  available, see \code{\link{disProg2sts}}.

  The constructor function \code{sts} can be used to setup an
  \code{"sts"} object.
  For areal time series, it can also capture a \code{map}
  of the regions, where the counts originate from.
  See Section \dQuote{Slots} below for a description of all class
  components, and Section \dQuote{Methods} for a list of extraction,
  conversion and visualization methods.

}

\usage{
sts(observed, start = c(2000, 1), frequency = 52,
    epoch = NULL, population = NULL, map = NULL, ...)
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
    of 2000. The \code{epoch} and \code{epochInYear} methods use the ISO
    8601 specification when converting between week numbers and dates,
    see \code{\link{isoWeekYear}}.}
  \item{epoch}{observation times, either as an integer sequence (default)
    or as a \code{Date} vector (in which case \code{epochAsDate} is
    automatically set to \code{TRUE}).}
  \item{population}{the population sizes (or fractions) underlying the
    \code{observed} counts, especially relevant for multivariate time
    series (for incidence maps or as offsets in epidemic models).
    If assumed constant over time, \code{population} can be supplied as
    a vector of length the number of columns (regions, groups, etc.) in
    \code{observed}. Otherwise, \code{population} needs to be a matrix
    of the same \code{\link{dim}}ension as \code{observed}. The
    \code{population} argument is an alias for the corresponding slot
    \code{populationFrac}. The default \code{NULL} value sets
    equal population fractions across all units.}
  \item{map}{optional spatial data representing the regions, either of
    class \code{"\linkSPclass{SpatialPolygons}"}
    (or \code{"\linkSPclass{SpatialPolygonsDataFrame}"}) or of class
    \code{"\link[sf]{sf}"} (requires package \CRANpkg{sf}).
    Its \code{row.names()} must contain the region IDs to be matched
    against \code{colnames(observed)}.}
  \item{\dots}{further named arguments with names corresponding to slot
    names (see the list below). For instance, in the public health surveillance context,
    the \code{state} slot is used to indicate outbreaks (default:
    \code{FALSE} for all observations). For areal time series data, the
    \code{map} and \code{neighbourhood} slots are used to store the
    spatial structure of the observation region.}
}

\section{Slots}{
  \describe{
    \item{\code{epoch}:}{a numeric vector specifying
      the time of observation, typically a week index. Depending on
      the \code{freq} slot, it could also index days or months.
      Furthermore, if \code{epochAsDate=TRUE} then \code{epoch}
      is the integer representation of \code{\link{Date}}s
      giving the exact date of the observation.}
    \item{\code{freq}:}{number of observations per year, e.g.,
      52 for weekly data, 12 for monthly data.}
    \item{\code{start}:}{vector of length two denoting the year and the
      sample number (week, month, etc.) of the first observation.}
    \item{\code{observed}:}{matrix of size \code{length(epoch)} times the
      number of regions containing the weekly/monthly number of counts in
      each region. The colnames of the matrix should match the ID values of
      the shapes in the \code{map} slot.}
    \item{\code{state}:}{matrix with the same dimensions as \code{observed}
      containing Booleans whether at the specific time point there was an
      outbreak in the region.}
    \item{\code{alarm}:}{matrix with the same dimensions as
      \code{observed} specifying whether an outbreak detection algorithm
      declared a specific time point in the region as having an alarm.}
    \item{\code{upperbound}:}{matrix with upper-bound values.}
    \item{\code{neighbourhood}:}{symmetric matrix of size
      \eqn{(number of regions)^2} describing the neighbourhood structure. It
      may either be a binary adjacency matrix or contain neighbourhood orders
      (see the Examples for how to infer the latter from the \code{map}).}
    \item{\code{populationFrac}:}{\code{matrix} of population 
      fractions or absolute numbers (see \code{multinomialTS} below)
      with dimensions \code{dim(observed)}.}
    \item{\code{map}:}{object of class \code{"\linkSPclass{SpatialPolygons}"}
      (or \code{"\linkSPclass{SpatialPolygonsDataFrame}"})
      providing a shape of the areas which are monitored or modelled.}
    \item{\code{control}:}{\code{\link{list}} of settings; this is a
      rather free data type to be returned by the surveillance algorithms.}
    \item{\code{epochAsDate}:}{a Boolean indicating
      if the \code{epoch} slot corresponds to \code{Date}s.}
    \item{\code{multinomialTS}:}{a Boolean
      stating whether to interpret the object as \code{observed} out of
      \code{population}, i.e. a multinomial interpretation instead of a
      count interpretation.}
  }
}

\section{Methods}{

  \subsection{Extraction of slots}{
  There is an extraction (and replacement) method for almost every slot.
  The name of the method corresponds to the slot name, with three exceptions:
  the \code{freq} slot can be extracted by \code{frequency()},
  the \code{populationFrac} slot is accessed by \code{population()},
  and the \code{alarm} slot is accessed by \code{alarms()}.

  \describe{
    \item{epoch}{\code{signature(x = "sts")}:
      extract the \code{epoch} slot. If the \code{sts} object is indexed
      by dates (\code{epochAsDate} = TRUE), the returned vector is of
      class \code{Date}, otherwise numeric (usually the integer
      sequence \code{1:nrow(x)}).\cr
      By explicitly requesting \code{epoch(x, as.Date = TRUE)}, dates
      can also be extracted if the \code{sts} object is not internally
      indexed by dates but has a standard frequency of 12 (monthly) or
      52 (weekly). The transformation is based on \code{start} and
      \code{freq} and will return the first day of each month
      (\code{freq=12}) and the Monday of each week (\code{freq=52}),
      respectively.}
    \item{frequency}{\code{signature(x = "sts")}:
      extract the \code{freq} slot.}
    \item{start}{\code{signature(x = "sts")}:
      extract the \code{start} slot.}
    \item{observed}{\code{signature(x = "sts")}:
      extract the \code{observed} slot.}
    \item{alarms}{\code{signature(x = "sts")}:
      extract the \code{alarm} slot.}
    \item{upperbound}{\code{signature(x = "sts")}:
      extract the \code{upperbound} slot.}
    \item{neighbourhood}{\code{signature(x = "sts")}:
      extract the \code{neighbourhood} slot.}
    \item{population}{\code{signature(x = "sts")}:
      extract the \code{populationFrac} slot.}
    \item{control}{\code{signature(x = "sts")}:
      extract the \code{control} slot.}
    \item{multinomialTS}{\code{signature(x = "sts")}:
      extract the \code{multinomialTS} slot.}
  }
  }

  \subsection{Other extraction methods}{
  \describe{
    \item{dim}{\code{signature(x = "sts")}:
      extract matrix dimensions of \code{observed}.
      This method also enables \code{nrow(x)} and \code{ncol(x)}.}
    \item{dimnames}{\code{signature(x = "sts")}:
      extract the \code{\link{dimnames}} of the \code{observed} matrix.
      This method also enables \code{rownames(x)} and \code{colnames(x)}.}
    \item{year}{\code{signature(x = "sts")}:
      extract the corresponding year of each observation.}
    \item{epochInYear}{\code{signature(x = "sts")}:
      extract the epoch number within the year.
      This corresponds to \code{\link[stats]{cycle}(as.ts(x))}.}

    \item{[}{\code{signature(x = "sts")}:
      subset rows (time points) and/or columns (units),
      see \code{help("\link{[,sts-method}")}.}
  }
  }

  \subsection{Transformation methods}{
  \describe{
    \item{aggregate}{\code{signature(x = "sts")}:
      see \code{\link{aggregate.sts}}.}

    \item{as.data.frame}{\code{signature(x = "sts")}:
      the default \code{as.data.frame} call will collect the following
      slots into a data frame: \code{observed}, \code{epoch},
      \code{state}, \code{alarm}, \code{upperbound}, and
      \code{populationFrac}. Additional columns will be created for
      \code{freq} (potentially varying by year for weekly or daily data
      if \code{x@epochAsDate} is \code{TRUE}) and
      \code{epochInPeriod} (the epoch fraction within the current year).\cr
      Calling the \code{as.data.frame} method with the argument
      \code{tidy = TRUE} will return \code{\link{tidy.sts}(x)},
      which reshapes multivariate \code{sts} objects to the
      \dQuote{long} format (one row per epoch and observational unit).
      The tidy format is particularly useful for standard regression
      models and customized plotting.}

    \item{coerce}{\code{signature(from="sts", to="ts")} and
      \code{signature(from="ts", to="sts")},
      to be called via \code{as(stsObj, "ts")} (or \code{as.ts(stsObj)})
      and \code{as(tsObj, "sts")}, respectively.}

    \item{as.xts}{convert to the \CRANpkg{xts} package format.}
  }
  }

  \subsection{Visualization methods}{
  \describe{
    \item{plot}{\code{signature(x = "sts", y = "missing")}:
      entry point to a collection of plot variants.
      The \code{type} of plot is specified using a formula,
      see \code{\link{plot.sts}} for details.}
    \item{autoplot}{a \CRANpkg{ggplot2} variant of the standard
      time-series-type plot, see \code{\link{autoplot.sts}}.}
    \item{animate}{see \code{\link{animate.sts}}.}
    \item{toLatex}{see \code{\link{toLatex.sts}}.}
  }
  }

}

\author{Michael \enc{Höhle}{Hoehle} and Sebastian Meyer}

\examples{
showClass("sts")

## create an sts object from time-series data
salmonellaDF <- read.table(system.file("extdata/salmonella.agona.txt",
                                       package = "surveillance"), header = TRUE)
str(salmonellaDF)
salmonella <- with(salmonellaDF,
                   sts(observed = observed, state = state,
                       start = c(1990, 1), frequency = 52))
salmonella
plot(salmonella)

## these data are also available as a legacy "disProg" object in the package
data(salmonella.agona)
stopifnot(all.equal(salmonella, disProg2sts(salmonella.agona)))


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
                        nbOrder(poly2adjmat(map))))
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
\keyword{ts}
\keyword{methods}
