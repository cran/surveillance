\name{polyAtBorder}
\alias{polyAtBorder}
\title{Indicate Polygons at the Border}

\description{
Determines which polygons of a \code{"\linkS4class{SpatialPolygons}"}
object are at the border, i.e. have coordinates in common with the
spatial union of all polygons (constructed using
\code{\link{unionSpatialPolygons}}).
}

\usage{
polyAtBorder(SpP, snap = sqrt(.Machine$double.eps),
             method = "rgeos", ...)
}

\arguments{
  \item{SpP}{
    an object of class \code{"\linkS4class{SpatialPolygons}"}.
  }
  \item{snap}{
    tolerance used to consider coordinates as identical.
  }
  \item{method}{method to use for \code{\link{unionSpatialPolygons}}.
    Defaults to \code{"rgeos"}, since \pkg{polyclip} uses integer
    arithmetic, which causes rounding errors usually requiring tuning of
    (i.e., increasing) the tolerance parameter \code{snap} (see example
    below).}
  \item{\dots}{further arguments passed to the chosen \code{method}.}
}

\value{
logical vector of the same length as \code{SpP} also inheriting its
\code{row.names}.
}

\author{
  Sebastian Meyer
}

\examples{
## Load districts of Germany
load(system.file("shapes", "districtsD.RData", package = "surveillance"))

## Determine districts at the border and check the result on the map
if (requireNamespace("rgeos")) {
    atBorder <- polyAtBorder(districtsD, method = "rgeos")
    plot(districtsD, col = atBorder)
}

## For method = "polyclip", a higher snapping tolerance is required
## to obtain the correct result
if (requireNamespace("polyclip")) {
    atBorder <- polyAtBorder(districtsD, snap = 1e-6, method = "polyclip")
    plot(districtsD, col = atBorder)
}
}

\keyword{spatial}
