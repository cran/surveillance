\name{unionSpatialPolygons}
\alias{unionSpatialPolygons}
\title{
  Compute the Unary Union of \code{"SpatialPolygons"}
}

\description{
Union all subpolygons of a
\code{"\link[sp:SpatialPolygons-class]{SpatialPolygons}"} object.
This is a wrapper for the polygon clipping engines implemented by
packages \pkg{rgeos} and \pkg{polyclip}.
}

\usage{
unionSpatialPolygons(SpP, method = c("rgeos", "polyclip"), ...)
}

\arguments{
  \item{SpP}{
    an object of class
    \code{"\link[sp:SpatialPolygons-class]{SpatialPolygons}"}.
    For the \pkg{polyclip} \code{method} only, all polygon classes for
    which an \code{\link{xylist}}-method exists should work as input.
  }
  \item{method}{
    polygon clipping machinery to use. Default is to simply call
    \code{\link[rgeos]{gUnaryUnion}} in package \pkg{rgeos}.
    For \code{method="polyclip"}, function
    \code{\link[polyclip]{polyclip}} from package \pkg{polyclip} is used.
    The alternative \code{method="gpclib"} via package \pkg{maptools}
    is no longer provided in \pkg{surveillance} > 1.20.3.
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
  \code{\link[rgeos]{gUnaryUnion}} in package \pkg{rgeos},
  \code{\link[polyclip]{polyclip}} in package \pkg{polyclip}.
}

\examples{
## Load districts of Germany
load(system.file("shapes", "districtsD.RData", package = "surveillance"))
plot(districtsD, border = "gray", asp = 1)  % to avoid sp -> sf

## Union these districts using either "rgeos" or "polyclip"
if (requireNamespace("rgeos"))  {
    stateD <- unionSpatialPolygons(districtsD, method = "rgeos")
    plot(stateD, add = TRUE, border = 2, lwd = 2)
}
if (requireNamespace("polyclip")) {
    stateD_pc <- unionSpatialPolygons(districtsD, method = "polyclip")
    plot(stateD_pc, add = TRUE, border = 1, lwd = 2, lty = 2)
}
}

\keyword{spatial}
