\name{intersectPolyCircle}
\alias{intersectPolyCircle}
\alias{intersectPolyCircle.owin}
\alias{intersectPolyCircle.SpatialPolygons}
\alias{intersectPolyCircle.gpc.poly}

\title{
  Intersection of a Polygonal and a Circular Domain
}

\description{
  This is a unifying wrapper around functionality of various packages
  dealing with spatial data. It computes the intersection of a circular
  domain and a polygonal domain (whose class defines the specific method).
}

\usage{
intersectPolyCircle(object, center, radius, ...)

\method{intersectPolyCircle}{owin}(object, center, radius, npoly = 32, ...)
\method{intersectPolyCircle}{SpatialPolygons}(object, center, radius, npoly = 32, ...)
\method{intersectPolyCircle}{gpc.poly}(object, center, radius, npoly = 32,
                   useGEOS = FALSE, ...)
}

\arguments{
  \item{object}{a polygonal domain of one of the supported classes.}
  \item{center,radius,npoly}{see \code{\link{discpoly}}.}
  \item{useGEOS}{logical indicating if package \pkg{rgeos}
    (\code{\link[rgeos]{gIntersection}}) should be used instead of
    package \pkg{gpclib}. The latter (default) requires explicit
    acceptance of \pkg{gpclib}'s restricted license via
    \code{\link{surveillance.options}(gpclib=TRUE)}.}
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
library("spatstat")
plot(letterR)
plot(intersectPolyCircle(letterR, c(3,2), 1), add=TRUE, col=2, lwd=3)
}

\keyword{spatial}
\keyword{manip}
