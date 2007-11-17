\name{sts-class}
\docType{class}
\alias{sts-class}
%
% other names
%\alias{aggregate,sts-method}
%\alias{plot,sts-method}
%\alias{plot,sts,missing-method}
\alias{nrow,sts-method}
\alias{ncol,sts-method}
\alias{dim,sts-method}
\alias{colnames,sts,missing,missing-method}
\alias{initialize,sts-method}
\encoding{latin1}

\title{Class "sts" -- surveillance time series}

\description{

  This is a rather leightweight class to implement multivariate time
  series of count used for public health surveillance data. The class
  captures the time series data as well as the spatial layout of the
  regions, where the data originate from.

}

\section{Slots}{
  \describe{
  \item{\code{week}:}{Object of class \code{numeric} specifying the week
  numbers. Actually this is not really used at the moment. }
  \item{\code{freq}:}{If weekly data \code{freq} corresponds to 52, in case of monthly data \code{freq} is 12.}
  \item{\code{start}:}{vector of length two denoting the year and the sample number (week, month, etc.) of the first observation}
  \item{\code{observed}:}{A matrix of size \code{length(week)} times the
  number of regions containing the weekly/monthly number of counts in
  each region. The colnames of the matrix should match the ID values of
  the shapes in the \code{map} slot.}
  \item{\code{state}:}{Matrix with the same dimension as \code{observed} containing booleans whether at the specific time point there was an outbreak in the region}
  \item{\code{alarm}:}{Matrix with the same dimension as \code{observed} specifying whether an outbreak detection algorithm declared a specific time point in the region as having an alarm. If the object containins just observations then this slot is null. }
  \item{\code{upperbound}:}{Matrix with upper bound values }
  \item{\code{neighbourhood}:}{Symmetric matrix of booleans size \eqn{(number of regions)^2} stating the neighbourhood matrix. }
  \item{\code{populationFrac}:}{Object of class \code{matrix}.}
  \item{\code{map}:}{Object of class \code{SpatialPolygonsDataFrame} providing a shape of the areas which are monitored. }
  \item{\code{control}:}{Object of class \code{list}, thais is a rather free data type to be returned by the surveillance algorithms. }
  }
}

\section{Methods}{
  \describe{ \item{nrow}{\code{signature(x = "sts")}:
      extract number of rows of the \code{observed} matrix slot. The
      dimension of the other matrix slots is similar.}
    \item{ncol}{\code{signature(x = "sts")}:
      extract number of columns of the \code{observed} matrix slot.}
    \item{dim}{\code{signature(x = "sts")}: extract matrix dimensions of
      \code{observed} using \code{\link{dim}}.}
    \item{colnames}{\code{signature(x="sts",do.NULL="missing",prefix="missing")}:
      extract \code{\link{colnames}} of the \code{observed} matrix.}
    \item{initialize}{\code{signature(x="sts")}: the internal function
      \code{init.sts} is called, which assigns all slots.  }
    \item{aggregate}{\code{signature(x="sts")}: see
      \code{\link{aggregate,sts-method}}
    } 
    \item{plot}{\code{signature(x="sts",y="missing",function(x, y, type,...)...)}: this function is
      the successor of the \code{\link{plot.disProg}} and 
      \code{\link{plot.survRes}} functions. It takes (more or less) the
      same arguments as \code{plot.survRes}. The most important
      difference is the type of plot, which is specified using
      \code{type}. See \code{\link{show,sts-method}} for details.
    }
  }
}

\author{M. Höhle}
\examples{
  data("ha")
  shp <- system.file("shapes/berlin.shp",package="surveillance")
  ha <- disProg2sts(ha, map=readShapePoly(shp,IDvar="SNAME"))
  plot(ha,type=observed ~ 1 | unit)

}
\keyword{classes}
