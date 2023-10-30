\name{intersectPolyCircle}
\alias{intersectPolyCircle}
\alias{intersectPolyCircle.owin}
%\alias{intersectPolyCircle.SpatialPolygons} # defunct in surveillance 1.22.0

\title{
  Intersection of a Polygonal and a Circular Domain
}

\description{
  This is a unifying wrapper around functionality of various packages
  dealing with spatial data. It computes the intersection of a circular
  domain and a polygonal domain (whose class defines the specific method).

  Currently the only supported class is \code{"\link[spatstat.geom]{owin}"}
  from package \CRANpkg{spatstat.geom}.
}

\usage{
intersectPolyCircle(object, center, radius, ...)

\method{intersectPolyCircle}{owin}(object, center, radius, npoly = 32, ...)
%\method{intersectPolyCircle}{SpatialPolygons}(object, center, radius, npoly = 32, ...)
}

\arguments{
  \item{object}{a polygonal domain of one of the supported classes.}
  \item{center,radius,npoly}{see \code{\link{discpoly}}.}
  \item{\dots}{potential further arguments (from the generic).}
}

\value{
  a polygonal domain of the same class as the input \code{object}.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{discpoly}} to generate a polygonal approximation to a disc
}

\examples{
library("spatstat.geom")
plot(letterR)
plot(intersectPolyCircle(letterR, c(3,2), 1), add=TRUE, col=2, lwd=3)
}

\keyword{spatial}
\keyword{manip}
