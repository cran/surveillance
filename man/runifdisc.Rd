\name{runifdisc}
\alias{runifdisc}

\title{
Sample Points Uniformly on a Disc
}

\description{
Sample \code{n} points uniformly on a disc of radius \code{r} in
two-dimensional euclidean space via transformation to polar coordinates:
the angle is sampled uniformly from \eqn{U(0,2\pi)}, the length is
sampled uniformly from \eqn{\sqrt{U(0,r^2)}}. The sampled polar
coordinates are then back-transformed to cartesian coordinates.
}

\usage{
runifdisc(n, r = 1)
}

\arguments{
  \item{n}{
    integer size of the sample.
  }
  \item{r}{
    numeric radius of the disc (centered at (0,0)).
  }
}

\value{
A two-column coordinate matrix of the sampled points.
}

\author{
Sebastian Meyer
}

\seealso{
\code{\link[spatstat]{runifdisc}} in package \pkg{spatstat}, which is
slightly more flexible and integrated within the \code{"ppp"} class.
}

\examples{
plot(surveillance:::runifdisc(1000, 3))
}

\keyword{datagen}
\keyword{distribution}
