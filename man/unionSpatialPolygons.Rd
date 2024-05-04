\name{unionSpatialPolygons}
\alias{unionSpatialPolygons}
\title{
  Compute the Unary Union of \code{"SpatialPolygons"}
}

\description{
Union all subpolygons of a
\code{"\link[sp:SpatialPolygons-class]{SpatialPolygons}"} object.
This is a legacy wrapper for the polygon clipping engines implemented by
packages \pkg{sf} and \pkg{polyclip}.
Internally, both \code{method}s need to convert the input polygons to a
class appropriate for the \code{method}, so are rather inefficient.
}

\usage{
unionSpatialPolygons(SpP, method = c("sf", "polyclip"), ...)
}

\arguments{
  \item{SpP}{
    an object of class
    \code{"\link[sp:SpatialPolygons-class]{SpatialPolygons}"}.
    For the \pkg{polyclip} \code{method} only, all polygon classes for
    which an \code{\link[polyCub]{xylist}}-method exists should work as input.
  }
  \item{method}{
    polygon clipping machinery to use. Default is to call
    \code{\link[sf]{st_union}} in package \pkg{sf}.
    For \code{method="polyclip"}, function
    \code{\link[polyclip]{polyclip}} from package \pkg{polyclip} is used.
    The old \code{method="gpclib"} is no longer available.
  }
  \item{\dots}{further arguments passed to the chosen \code{method}.}
}

\value{
  an object of class
  \code{"\link[sp:SpatialPolygons-class]{SpatialPolygons}"} representing
  the union of all subpolygons.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link[sf]{st_union}} in package \pkg{sf},
  \code{\link[polyclip]{polyclip}} in package \pkg{polyclip}.
}

\examples{
## Load districts of Germany
load(system.file("shapes", "districtsD.RData", package = "surveillance"))
plot(districtsD, border = "gray", asp = 1)  % to avoid sp -> sf

## Union these districts using either "sf" or "polyclip"
if (requireNamespace("sf"))  {
    stateD <- unionSpatialPolygons(districtsD, method = "sf")
    plot(stateD, add = TRUE, border = 2, lwd = 2)
}
if (requireNamespace("polyclip")) {
    stateD_pc <- unionSpatialPolygons(districtsD, method = "polyclip")
    plot(stateD_pc, add = TRUE, border = 1, lwd = 2, lty = 2)
}
}

\keyword{spatial}
