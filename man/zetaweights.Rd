\name{zetaweights}
\alias{zetaweights}
\title{
  Power-Law Weights According to Neighbourhood Order
}

\description{
  Compute power-law weights with decay parameter \code{d}
  based on a matrix of neighbourhood orders \code{nbmat}
  (e.g., as obtained via \code{\link{nbOrder}}).
  Without normalization and truncation,
  this is just \eqn{o^{-d}} (where \eqn{o} is a neighbourhood order).
  This function is mainly used internally for \code{\link{W_powerlaw}}
  weights in \code{\link{hhh4}} models.
}

\usage{
zetaweights(nbmat, d = 1, maxlag = max(nbmat), normalize = FALSE)
}

\arguments{
  \item{nbmat}{numeric, symmetric matrix of neighbourhood orders.}
  \item{d}{single numeric decay parameter (default: 1). Should be
  positive.}
  \item{maxlag}{single numeric specifying an upper limit for the power
    law. For neighbourhood orders > \code{maxlag}, the resulting weight
    is 0. Defaults to no truncation.}
  \item{normalize}{Should the resulting weight matrix be normalized such
    that rows sum to 1?}
}

\value{
  a numeric matrix with same dimensions and names as the input matrix.
}

\author{
  Sebastian Meyer
}

\seealso{\code{\link{W_powerlaw}}}

\examples{
nbmat <- matrix(c(0,1,2,2,
                  1,0,1,1,
                  2,1,0,2,
                  2,1,2,0), 4, 4, byrow=TRUE)
zetaweights(nbmat, d=1, normalize=FALSE) # harmonic: o^-1
zetaweights(nbmat, d=1, normalize=TRUE)  # rowSums=1
zetaweights(nbmat, maxlag=1, normalize=FALSE) # results in adjacency matrix
}

\keyword{spatial}
\keyword{utilities}
