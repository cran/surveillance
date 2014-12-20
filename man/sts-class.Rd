\name{sts-class}
\Rdversion{1.1}
\docType{class}
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

\encoding{latin1}

\title{Class \code{"sts"} -- surveillance time series}

\description{

  This is a lightweight S4 class to implement multivariate time
  series of counts used for public health surveillance data. The class
  captures the time series data as well as the spatial layout of the
  regions, where the data originate from.

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
      declared a specific time point in the region as having an alarm. If
      the object containins just observations then this slot is null. }
    \item{\code{upperbound}:}{Matrix with upper bound values }
    \item{\code{neighbourhood}:}{Symmetric matrix of size
      \eqn{(number of regions)^2} describing the neighbourhood structure. It
      may either be a binary adjacency matrix or contain neighbourhood orders.}
    \item{\code{populationFrac}:}{A \code{matrix} of population 
      fractions (with dimensions \code{dim(observed)}).}
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
      \code{observed} using \code{\link{dim}}.}
    \item{observed}{\code{signature(x = "sts")}: extract the \code{observed}
      slot of an \code{sts} object.}
    \item{population}{\code{signature(x = "sts")}: extract the \code{population}
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
    \item{dimnames}{\code{signature(x="sts")}: extract the \code{\link{dimnames}} of the \code{observed} matrix.}
    \item{initialize}{\code{signature(x="sts")}: the internal function
      \code{init.sts} is called, which assigns all slots.  }
    \item{aggregate}{\code{signature(x="sts")}: see
      \code{\link{aggregate,sts-method}}}
    \item{year}{\code{signature(x = "sts")}:
      extracts the corresponding year of each observation of \code{x}}
    \item{obsinyear}{\code{signature(x = "sts")}:
      extracts the corresponding week number within the year of each
      observation of \code{x}}
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
  }
}

\author{M. \enc{Höhle}{Hoehle}}
\examples{
if (requireNamespace("maptools")) {
   # load disProg-object "ha" and convert to S4-class "sts"
   data("ha")
   shpfile <- system.file("shapes/berlin.shp",package="surveillance")
   ha.sts <- disProg2sts(ha, map=maptools::readShapePoly(shpfile,IDvar="SNAME"))
} else {
   data("ha.sts")
   # is almost identical to the above except that German umlauts
   # have been replaced in 'ha.sts@map@data$BEZIRK' for compatibility reasons
}

ha.sts
plot(ha.sts, type=observed ~ 1 | unit)

#Convert ts/mts object to sts
z <- ts(matrix(rpois(300,10), 100, 3), start = c(1961, 1), frequency = 12)
sts.z <- as(z,"sts")
plot(sts.z)

}

\keyword{classes}
