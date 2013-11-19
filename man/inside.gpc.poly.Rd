\name{inside.gpc.poly}
\alias{inside.gpc.poly}
\title{
Test Whether Points are Inside a \code{"gpc.poly"} Polygon
}

\description{
Same as, e.g., \code{\link[spatstat]{inside.owin}} from package \pkg{spatstat}
and \code{\link[sp]{point.in.polygon}} from package \pkg{sp},
i.e., test whether points lie inside or outside a given polygon.
Actually, the method for \code{"gpc.poly"} documented here internally uses
the \code{\link[sp]{point.in.polygon}} function.
}

\usage{
inside.gpc.poly(x, y = NULL, polyregion, mode.checked = FALSE)
}

\arguments{
  \item{x,y}{
    numeric vectors of coordinates of the points to be tested.
    The coordinates can be supplied in any form accepted by
    \code{\link{xy.coords}}.
}
  \item{polyregion}{
    an object of class \code{"gpc.poly"}. It is checked if the points specified
    through \code{x} and \code{y} fall into this polygonal region.
}
  \item{mode.checked}{
    passed to \code{\link[sp]{point.in.polygon}}.
}
}

\details{
The nodes and edges of (non-hole) polygons are treated as being inside.
Points that fall \emph{strictly} inside holes are treated as being outside
of the polygon.
}

\value{
Logical vector whose \code{i}th entry is \code{TRUE} if the corresponding
point \code{(x[i],y[i])} is inside \code{polyregion}.
}

\author{
Sebastian Meyer
}

\examples{
if (requireNamespace("rgeos")) {
    poly <- discpoly(c(0.5,0.5), 0.5, npoly=4, class="gpc.poly")
    pts <- cbind(x=runif(50), y=runif(50))
    plot(poly)
    points(pts, col=1+inside.gpc.poly(pts, polyregion=poly))
}
}

\keyword{utilities}
\keyword{spatial}
