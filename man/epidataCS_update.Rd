\name{epidataCS_update}
\alias{update.epidataCS}
\title{
  Update method for \code{"epidataCS"}
}

\description{
The \code{\link{update}} method for the \code{"\link{epidataCS}"} class
may be used to modify the hyperparameters \eqn{\epsilon} (\code{eps.t})
and \eqn{\delta} (\code{eps.s}), the indicator matrix \code{qmatrix} determining
possible transmission between the event types, the numerical
accuracy \code{nCircle2Poly} of the polygonal approximation, and
the endemic covariates from \code{stgrid} (including the time intervals).
The update method will also update the auxiliary information contained
in an \code{"epidataCS"} object accordingly, e.g., the vector of potential
sources of each event, the influence regions, or the endemic covariates
copied from the new \code{stgrid}.
}

\usage{
\method{update}{epidataCS}(object, eps.t, eps.s, qmatrix, nCircle2Poly, stgrid, ...)
}

\arguments{
  \item{object}{
    an object of class \code{"epidataCS"}.
}
  \item{eps.t}{
    numeric vector of length 1 or corresponding to the number of events in
    \code{object$events}.  The event data column \code{eps.t} specifies
    the maximum temporal influence radius (e.g., length of infectious
    period, time to culling, etc.) of the events. 
}
  \item{eps.s}{
    numeric vector of length 1 or corresponding to the number of events in
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
  \item{stgrid}{
    a new \code{data.frame} with endemic covariates, possibly
    transformed from or adding to the original \code{object$stgrid}.
    The grid must cover the same regions as the original, i.e.,
    \code{levels(object$stgrid$tile)} must remain identical.
    See \code{\link{epidataCS}} for a detailed description of the
    required format.
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

## use a time-constant grid
imdepi3 <- update(imdepi, stgrid = subset(imdepi$stgrid, BLOCK == 1, -stop))
(s3 <- summary(imdepi3)) # "1 time block"
}

\keyword{manip}
\keyword{utilities}
\keyword{methods}
