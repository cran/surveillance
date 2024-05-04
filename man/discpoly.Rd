\name{discpoly}
\alias{discpoly}

\title{Polygonal Approximation of a Disc/Circle}

\description{
  Generates a polygon representing a disc/circle (in planar
  coordinates) as an object of one of three possible
  classes: \code{"\link[sp:Polygon-class]{Polygon}"} from package \CRANpkg{sp},
  \code{"\link[spatstat.geom]{owin}"} from package \CRANpkg{spatstat.geom}, or
  \code{"gpc.poly"} from \pkg{gpclib} (if available).
}

\usage{
discpoly(center, radius, npoly = 64,
         class = c("Polygon", "owin", "gpc.poly"),
         hole = FALSE)
}

\arguments{
  \item{center}{numeric vector of length 2 (center coordinates of the circle).}

  \item{radius}{single numeric value (radius of the circle).}

  \item{npoly}{single integer. Number of edges of the polygonal approximation.}

  \item{class}{class of the resulting polygon (partial name
  matching applies). For \code{"owin"}, this is just a
  wrapper around \pkg{spatstat.geom}'s own \code{\link[spatstat.geom]{disc}} function.}

  \item{hole}{logical. Does the resulting polygon represent a hole?}
}

\value{
  A polygon of class \code{class} representing a
  circle/disc with \code{npoly} edges accuracy.
  
  If \code{class="gpc.poly"} and this S4 class is not yet registered
  in the current \R session (by loading \pkg{gpclib} beforehand), only the
  \code{pts} slot of a \code{"gpc.poly"} is returned with a warning.
}

\examples{
## Construct circles with increasing accuracy and of different spatial classes
disc1 <- discpoly(c(0,0), 5, npoly=4, class = "owin")
disc2 <- discpoly(c(0,0), 5, npoly=16, class = "Polygon")
disc3 <- discpoly(c(0,0), 5, npoly=64, class = "gpc.poly")  # may warn

## Look at the results
print(disc1)
plot(disc1, axes=TRUE, main="", border=2)

str(disc2)
lines(disc2, col=3)

str(disc3)  # a list or a formal "gpc.poly" (if gpclib is available)
if (is(disc3, "gpc.poly")) {
  plot(disc3, add=TRUE, poly.args=list(border=4))
} else {
  lines(disc3[[1]], col=4)
}

## to only _draw_ a circle
symbols(0, 0, circles=5, inches=FALSE, add=TRUE, fg=5)
}

\seealso{
  \code{\link[spatstat.geom]{disc}} in package \pkg{spatstat.geom}.
}

\keyword{datagen}
\keyword{spatial}
