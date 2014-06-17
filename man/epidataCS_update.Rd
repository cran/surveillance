\name{epidataCS_update}
\alias{update.epidataCS}
\title{
  Update method for \code{"epidataCS"}
}

\description{
The \code{\link{update}} method for the \code{"\link{epidataCS}"} class
may be used to modify the hyperparameters \eqn{\epsilon} (\code{eps.t})
and \eqn{\delta} (\code{eps.s}), the indicator matrix \code{qmatrix} of
possible ways of transmission between the event types, and the numerical
accuracy \code{nCircle2Poly} of the polygonal representation of a
circle.
The update method will also update the auxiliary information contained
in an \code{"epidataCS"} object accordingly, e.g., the vector of potential
sources of each event, or the polygonal representation of the influence
region.
}

\usage{
\method{update}{epidataCS}(object, eps.t, eps.s, qmatrix, nCircle2Poly, ...)
}

\arguments{
  \item{object}{
    an object of class \code{"epidataCS"}.
}
  \item{eps.t}{
    numeric vector of length the number of events in
    \code{object$events}.  The event data column \code{eps.t} specifies
    the maximum temporal influence radius (e.g., length of infectious
    period, time to culling, etc.) of the events. 
}
  \item{eps.s}{
    numeric vector of length the number of events in
    \code{object$events}.  The event data column \code{eps.s} specifies
    the maximum spatial influence radius of the events.
}
  \item{qmatrix}{
    square indicator matrix (0/1 or TRUE/FALSE) for possible
    transmission between the event types.
}
  \item{nCircle2Poly}{
    accuracy (number of edges) of the polygonal approximation of a circle.
}
  \item{\dots}{
    unused (argument of the generic).
}
}

\value{
  The updated \code{"epidataCS"} object.
}

\author{
  Sebastian Meyer
}

\seealso{
  class \code{"\link{epidataCS}"}.
}

\examples{
data("imdepi")

## assume different interaction ranges and simplify polygons
imdepi2 <- update(imdepi, eps.t = 20, eps.s = Inf, nCircle2Poly = 16)
    
(s <- summary(imdepi))
(s2 <- summary(imdepi2))
## The update reduced the number of infectives (along time)
## because the length of the infectious periods is reduced. It also 
## changed the set of potential sources of transmission for each
## event, since the interaction is shorter in time but wider in space
## (eps.s=Inf means interaction over the whole observation region).
}

\keyword{manip}
\keyword{utilities}
\keyword{methods}
