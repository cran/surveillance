\encoding{latin1}
\name{polyCub}
\alias{polyCub}
\alias{polyCub.midpoint}
\alias{polyCub.SV}
\alias{polyCub.exact.Gauss}

\title{
  Two-Dimensional Numerical Integration over a Polygonal Domain
}

\description{
The functions documented here numerically integrate a two-dimensional
function over a polygonal domain. The first method uses the classical midpoint rule
converting the polygons to binary pixel images using the
\code{\link[spatstat]{as.im.function}} method from package \pkg{spatstat}.
The second one uses the more sophisticated product Gauss cubature
proposed by Sommariva & Vianello (2007).
The third one is a quasi-exact method specific to the integration of the
bivariate Gaussian density over polygonal domains and is based on formulae
from Chapter 26 of the famous Abramowitz & Stegun (1970) handbook, i.e.
triangulation of the polygonal domain and appropriate evaluations of
\code{\link[mvtnorm]{pmvnorm}} from package \pkg{mvtnorm}.
See Section 3.2 of Meyer (2010) for a more detailed description and benchmark
experiment of those three (among others).
}

\usage{
polyCub(polyregion, f, method = c("SV", "midpoint", "exact.Gauss"), ...,
        plot = FALSE)

polyCub.midpoint(polyregion, f, ...,
                 eps = NULL, dimyx = NULL, plot = FALSE)
polyCub.SV(polyregion, f, ...,
           nGQ = 20, alpha = NULL, rotation = FALSE, plot = FALSE)
polyCub.exact.Gauss(polyregion, mean = c(0,0), Sigma = diag(2),
                    plot = FALSE)
}

\arguments{
  \item{polyregion}{
    a polygon representing the integration domain.
    For \code{polyCub.midpoint}, \code{polyregion} can be anything coercible to the
    \pkg{spatstat} class \code{"\link[spatstat]{owin}"} (by
    \code{\link[spatstat]{as.owin}}).
    The \code{polyCub.SV} method allows for classes \code{"\link[spatstat]{owin}"},
    \code{"\link[gpclib:gpc.poly-class]{gpc.poly}"}, and
    \code{"\linkS4class{SpatialPolygons}"}.

    Finally, the \code{polyCub.exact.Gauss} method assumes
    \code{polyregion} to be anything coercible to a
    \code{"\link[gpclib:gpc.poly-class]{gpc.poly}"} polygon.\cr
    Unfortunately, it seems that there is no convention for the implementation of
    polygons of class \code{"gpc.poly"}. Here we assume vertices
    ordered according to the \pkg{sp} convention, i.e. \emph{clockwise} for normal
    boundaries and \emph{anticlockwise} for holes. However, in contrast to \pkg{sp},
    the first vertex should \emph{not} be repeated!
    The package \pkg{surveillance} offers methods to convert
    \code{"SpatialPolygons"} and \code{"owin"} objects to class
    \code{"gpc.poly"}, which follow this convention. However,
    \pkg{gpclib} has a restricted license, which has to be accepted
    explicitly via setting \code{\link{surveillance.options}(gpclib=TRUE)}
    prior to using \code{polyCub.exact.Gauss}.
}
  \item{f}{
    two-dimensional function to be integrated. As its first argument the function must take
    a coordinate matrix, i.e. a matrix with two columns. For the \code{polyCub.exact.Gauss}
    method, \code{f} is ignored since it is specific to the bivariate normal density.
}
  \item{method}{
    choose one of the implemented cubature methods (partial argument matching is applied).
    Defaults to using the product Gauss cubature implemented in \code{polyCub.SV}.
  }
  \item{\dots}{
    further arguments passed to \code{f}.
}
  \item{eps}{
    numeric scalar specifying the (approximate) width/height of the
    pixels for the generation of the image.
}
  \item{dimyx}{
    vector of length two specifying the number of subdivisions in each dimension for the
    generation of the image. See \code{\link{as.im.function}} and \code{\link{as.mask}}.
}
  \item{nGQ}{
    degree of the one-dimensional Gauss-Legendre quadrature rule
    (default: 20). See \code{\link[statmod]{gauss.quad}} in package
    \pkg{statmod}, on which this function depends.
}
  \item{alpha}{
    base-line of the (rotated) polygon at \eqn{x = alpha} (see Sommariva
    & Vianello (2007) for an explication).
    If \code{NULL} (default), the midpoint of the x-range of the polygon is chosen
    if no \code{rotation} is performed, and otherwise the \eqn{x}-coordinate
    of the rotated point \code{"P"}.
    If \code{f} has its maximum value at the origin (0,0), which is the
    case for the typical usage in this package (spatial interaction
    functions like \code{\link{siaf.gaussian}}), \code{alpha = 0} is a
    reasonable choice.
  }
  \item{rotation}{
    logical (default: \code{FALSE}) or a list of points \code{"P"} and
    \code{"Q"} describing the preferred direction. If \code{TRUE}, the
    polygon is rotated according to the vertices \code{"P"} and
    \code{"Q"}, which are farthest apart (see Sommariva & Vianello,
    2007). For convex polygons, this rotation guarantees that all nodes
    fall inside the polygon.
  }
  \item{mean, Sigma}{
    mean and covariance matrix of the bivariate normal density to be integrated.
  }
  \item{plot}{
    logical indicating if an illustrative plot of the numerical integration should be produced.
  }
}

\value{
  The approximated value of the integral of \code{f} over \code{polyregion}.
  
  For the \code{polyCub.exact.Gauss} method, two attributes are appended to the integral value:
  \describe{
    \item{nEval}{number of triangles over which the standard bivariate normal density had to be integrated,
      i.e. number of calls to \code{\link[mvtnorm]{pmvnorm}} and \code{\link[stats]{pnorm}}, the former of
      which being the most time-consuming operation.
    }
    \item{error}{
      Approximate absolute integration error steming from the error introduced by the \code{nEval}
      \code{\link[mvtnorm]{pmvnorm}} evaluations. For this reason, the cubature method is in fact only
      quasi-exact (as is the \code{pmvnorm} function).
    }
  }
}

\references{
  Abramowitz, M. and Stegun, I. A. (1970). Handbook of Mathematical Functions with Formulas,
  Graphs, and Mathematical Tables (9th ed.). New York: Dover Publications.

  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{ä}{ae}t
  M\enc{ü}{ue}nchen.\cr
  Available as \url{http://epub.ub.uni-muenchen.de/11703/}
  
  Sommariva, A. and Vianello, M. (2007): Product Gauss cubature over polygons
  based on Green's integration formula. Bit Numerical Mathematics, 47 (2), 441-453.
}

\author{
  Sebastian Meyer;
  the original MATLAB implementation on which the product Gauss cubature in
  \code{polyCub.SV} is based, is due to Sommariva & Vianello (2007)
}

\note{
  The \R implementation of the product Gauss cubature in \code{polycub.SV} is
  based on the corresponding MATLAB code (\code{polygauss}) by Sommariva &
  Vianello (2007); but it only covers the case of \dQuote{simple} polygons. The
  efficiency rotation also proposed by Sommariva & Vianello (2007) is not
  implemented here since it does in general only apply to convex polygons.
}

\examples{
# two-dimensional function to integrate (here: isotropic Gaussian density with zero mean)
f <- function (s, sigma = 5) exp(-rowSums(s^2)/2/sigma^2) / (2*pi*sigma^2)

# simple polygonal integration domain
disc.owin <- discpoly(c(3,2), 5, npoly=8, class="owin")

# plot image of the function and integration domain
gr <- seq(-8,8,by=0.2)
vals <- matrix(f(expand.grid(gr,gr)), nrow=length(gr))
image(x=gr, y=gr, z=vals)
library("spatstat")    
plot.owin(disc.owin, hatch=TRUE, add=TRUE)

# Two-dimensional midpoint rule
testmidpoint <- function (eps, main=paste("Two-dimensional midpoint rule with eps =",eps))
{
    image(x=gr, y=gr, z=vals, main=main)
    plot.owin(disc.owin, hatch=TRUE, add=TRUE)
    # add binary mask to plot
    #plot(as.mask(disc.owin, eps=eps), add=TRUE, box=FALSE)
    # add evaluation points to plot
    with(as.mask(disc.owin, eps=eps), points(expand.grid(xcol, yrow), col=m, pch=20))
    polyCub.midpoint(disc.owin, f, eps=eps)
}
testmidpoint(5)
testmidpoint(3)
testmidpoint(0.5)
testmidpoint(0.2)

# Product Gauss cubature by Sommariva & Vianello (2007)
if (require("statmod")) {
    for (nGQ in c(1:5,10,20,60)) {
        cat("nGQ =", sprintf("\%2i",nGQ), ": ",
            format(polyCub.SV(disc.owin, f, nGQ=nGQ), digits=16),
            "\n")
    }
}

# quasi-exact cubature by triangulation (gpclib) and using mvtnorm::pmvnorm()
if (require("mvtnorm") && require("gpclib")) {
    oopt <- surveillance.options(gpclib=TRUE)
    print(polyCub.exact.Gauss(disc.owin, mean=c(0,0), Sigma=5^2*diag(2),
                              plot=TRUE), digits=16)
    surveillance.options(oopt)
}
}

\keyword{math}
\keyword{utilities}
\keyword{spatial}

