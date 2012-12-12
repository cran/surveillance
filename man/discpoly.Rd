\name{discpoly}
\alias{discpoly}
\title{
  Generate a Polygon Representing a Disc/Circle (in Planar Coordinates)
}

\description{
  This functions generates a circle/disc represented as a polygon of one of three
  possible data classes: \code{"\link[gpclib:gpc.poly-class]{gpc.poly}"},
  \code{"\link[spatstat]{owin}"}, and \code{"\link[sp:Polygon-class]{Polygon}"}.
}

\usage{
discpoly(center, r, npoly = 64, class = c("Polygon", "owin", "gpc.poly"),
         hole = FALSE)
}

\arguments{
  \item{center}{
    numeric vector of length two. Center coordinates of the circle.
  }
  \item{r}{
    numeric scalar. Radius of the circle.
  }
  \item{npoly}{
    integer. Number of edges of the polygonal approximation to the circle.
  }
  \item{class}{
    Which class of polygon should be generated? For \pkg{spatstat}'s
    \code{"owin"} class, the \link[spatstat]{disc} function is used.
  }
  \item{hole}{
    logical. \code{hole} flag of the polygon.
  }
}

\value{
  Polygon of class \code{class} representing a circle/disc with
  \code{npoly} edges accuracy. 
}

\author{
  Sebastian Meyer
}

\seealso{
  \link[spatstat]{disc} in package \pkg{spatstat}.
}

\examples{
# Construct circle with increasing accuracy and of different spatial classes
disc1 <- discpoly(c(0,0), 5, npoly=4, class = "owin")
disc2 <- discpoly(c(0,0), 5, npoly=16, class = "Polygon")

library("spatstat")
print(disc1)
plot(disc1, axes=TRUE, main="", border=2)

print(disc2)
lines(disc2, col=3)

if (require("gpclib")) {
    oopt <- surveillance.options(gpclib=TRUE)
    disc3 <- discpoly(c(0,0), 5, npoly=64, class = "gpc.poly")
    print(disc3)
    plot(disc3, add=TRUE, poly.args=list(border=4))
    surveillance.options(oopt)
}

# if one only wants to draw a circle without an object behind
symbols(0, 0, circles=5, inches=FALSE, add=TRUE, fg=5)
}

\keyword{spatial}
\keyword{datagen}
