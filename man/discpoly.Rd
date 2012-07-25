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
discpoly(center, r, npoly = 64, class = c("gpc.poly", "owin", "Polygon"),
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
    Which class of polygon should be generated?
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
disc1 <- discpoly(c(0,0), 5, npoly=4, class = "gpc.poly")
disc2 <- discpoly(c(0,0), 5, npoly=16, class = "owin")
disc3 <- discpoly(c(0,0), 5, npoly=64, class = "Polygon")

print(disc1)
plot(disc1)

library("spatstat")
print(disc2)
plot(disc2, add=TRUE, border=2)

print(disc3)
lines(disc3, col=3)

# if one only wants to draw a circle one could just use
symbols(0,0, circles=5, inches=FALSE, add=TRUE, fg=4)
}

\keyword{spatial}
\keyword{datagen}
