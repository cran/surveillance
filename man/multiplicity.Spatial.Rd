\name{multiplicity.Spatial}
\alias{multiplicity.Spatial}

\title{
Count Number of Instances of Points
}

\description{
  The generic function \code{multiplicity} defined in \pkg{spatstat} is
  intended to count the number of duplicates of each element of an object.
  \pkg{spatstat} already offers methods for point patterns, matrices and
  data frames, and here we add a method for \code{Spatial} objects from
  the \pkg{sp} package. It is a wrapper for the default method, which
  effectively computes the distance matrix of the points,
  and then just counts the number of zeroes in each row.
}

\usage{
\method{multiplicity}{Spatial}(x)
}

\arguments{
  \item{x}{
    a \code{"\linkS4class{Spatial}"} object (we only need a
    \code{\link{coordinates}}-method), e.g. of class
    \code{"\linkS4class{SpatialPoints}"}. 
  }
}

\value{
  an integer vector containing the number of instances of each
  point of the object.
}

\seealso{
  \code{\link[spatstat]{multiplicity}} in package \pkg{spatstat}.
  See the Examples of the \code{\link{hagelloch}} data for a specific
  use of \code{multiplicity}.
}

\examples{
foo <- SpatialPoints(matrix(c(1,2,
                              2,3,
                              1,2,
                              4,5), 4, 2, byrow=TRUE))
multiplicity(foo)

# the following function determines the multiplicities in a matrix
# or data frame and returns unique rows with appended multiplicity
countunique <- function(x) unique(cbind(x, count=multiplicity(x)))
countunique(coordinates(foo))
}

\keyword{utilities}
\keyword{spatial}
