\name{scale.gpc.poly}
\alias{scale.gpc.poly}

\title{Centering and Scaling a \code{"gpc.poly"} Polygon}

\description{
This is a re-implementation of the corresponding method from package
\pkg{gpclib} to also allow centering.
}

\usage{
\method{scale}{gpc.poly}(x, center = c(0,0), scale = c(1,1))
}

\arguments{
  \item{x}{an object of class \code{"gpc.poly"}.}
  \item{center}{numeric vector of length 2 (x,y), which will be subtracted
    from the respective coordinates of \code{x}.}
  \item{scale}{numeric vector of length 2 (x,y), which serves as the
    divisor for the respective coordinates of \code{x}.}
}

\value{
A \code{"gpc.poly"}, the shifted and/or scaled version of \code{x}.
}

\keyword{methods}
\keyword{manip}
