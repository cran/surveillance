\name{untie}
\alias{untie}
\alias{untie.epidataCS}
\alias{untie.matrix}
\alias{untie.default}
\title{
  Randomly Break Ties in Data
}

\description{
  This is a generic function intended to randomly break tied data in a
  way similar to what \code{\link{jitter}} does: tie-breaking is
  performed by shifting \emph{all} data points by a random amount.
  The \pkg{surveillance} package defines methods for matrices,
  \code{"epidataCS"}, and a default method for numeric vectors.
}

\usage{
untie(x, amount, ...)

\method{untie}{epidataCS}(x, amount = list(t=NULL, s=NULL),
      minsep = list(t=0, s=0), direction = "left", keep.sources = FALSE,
      ..., verbose = FALSE)
\method{untie}{matrix}(x, amount = NULL, minsep = 0,
      constraint = NULL, giveup = 1000, ...)
\method{untie}{default}(x, amount = NULL, minsep = 0,
      direction = c("symmetric", "left", "right"), sort = NULL,
      giveup = 1000, ...)
}

\arguments{
  \item{x}{
    the data to be untied.
  }
  \item{amount}{
    upper bound for the random amount by which data are shifted.
    \code{NULL} means to use a data-driven default, which equals
    the minimum separation of the data points for the non-symmetric
    default method and its half for the symmetric default method and the
    \code{matrix} method. 
  }
  \item{minsep}{minimum separation of jittered points. Can only be
    obeyed if much smaller than \code{amount} (also depending on the
    number of points). \code{minsep>0} is currently only implemented for
    the spatial (matrix) method.}
  \item{keep.sources}{
    logical (\code{FALSE}). If \code{TRUE}, the original list of
    possible event sources in \code{x$events$.sources} will be
    preserved. For instance, events observed at the same time did by
    definition not trigger each other; however, after random
    tie-breaking one event will precede the other and considered as a
    potential source of infection for the latter, although it could just
    as well be the other way round. Enabling \code{keep.sources} will
    use the \code{.sources} list from the original (tied)
    \code{"epidataCS"} object. Note, however, that an update is forced
    within \code{twinstim} if a subset of the data is selected
    for model fitting or if a different \code{qmatrix} is supplied.
  }
  \item{constraint}{
    an object of class \code{"\linkS4class{SpatialPolygons}"}
    representing the domain which the points of the matrix should belong
    to -- before and after jittering.
  }
  \item{giveup}{number of attempts after which the
    algorithm should stop trying to generate new points.}
  \item{direction}{
    one of \code{"symmetric"} (default), \code{"left"}, or
    \code{"right"}, indicating in which direction vector elements should
    be shifted.
  }
  \item{sort}{
    logical indicating if the jittered vector should be sorted. Defaults
    to doing so if the original vector was already sorted.
  }
  \item{\dots}{
    For the \code{"epidataCS"}-method: arguments passed to the
    \code{matrix}- or \code{default}-method (\code{giveup}).
    Unused in other methods.
  }
  \item{verbose}{logical passed to \code{\link{as.epidataCS}}.}
}

\details{
  For numeric vectors (default method), the jittered version is the
  same as for \code{\link{jitter}(x, amount=amount)}, if
  \code{direction="symmetric"} (and \code{amount} is non-\code{NULL}),
  and otherwise uses
  \code{x} \dQuote{+-} \code{runif(length(x), 0, amount)}.
  
  For matrices, a vector uniformly drawn from the disc with radius
  \code{amount} is added to each point (row).
  
  For \code{"epidataCS"}, \code{amount} is a list stating the amounts
  for the temporal and/or spatial dimension, respectively. It then
  uses the specific methods with arguments \code{constraint=x$W},
  \code{direction}, and \code{sort=TRUE}. Note that this implements a
  simplistic approach of tie-breaking where all events are assumed to be
  subject to the same amounts of censoring, and the default amounts may
  not be sensible choices.
}

\value{
  the untied (jittered) data.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{jitter}}
}

\examples{
# vector example
set.seed(123)
untie(c(rep(1,3), rep(1.2, 4), rep(3,3)), direction="left", sort=FALSE)

# spatial example
data(imdepi)
coords <- coordinates(imdepi$events)
table(duplicated(coords))
plot(coords, cex=sqrt(multiplicity(coords)))
set.seed(1)
coords_untied <- untie(coords)
stopifnot(!anyDuplicated(coords_untied))
points(coords_untied, col=2) # shifted by very small amount in this case
}

\keyword{utilities}
\keyword{manip}
\keyword{dplot}
