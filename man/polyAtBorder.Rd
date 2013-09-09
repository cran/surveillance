\name{polyAtBorder}
\alias{polyAtBorder}
\title{Indicate Polygons at the Border}

\description{
Determines which polygons of a \code{"\linkS4class{SpatialPolygons}"}
object are at the border, i.e. have coordinates in common with the
spatial union of all polygons (constructed using the function
\code{\link[maptools]{unionSpatialPolygons}} from package \pkg{maptools}).
}

\usage{
polyAtBorder(SpP, snap = sqrt(.Machine$double.eps))
}

\arguments{
  \item{SpP}{
    an object of class \code{"\linkS4class{SpatialPolygons}"}.
  }
  \item{snap}{
    tolerance used to consider coordinates as identical.
  }
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
load(system.file("shapes", "districtsD.RData", package="surveillance"))

## indicate districts at the border
districtAtBorder <- polyAtBorder(districtsD)

## plot to check
plot(districtsD, col=districtAtBorder)
}

\keyword{spatial}
