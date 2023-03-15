\name{nbOrder}
\alias{nbOrder}
\title{
  Determine Neighbourhood Order Matrix from Binary Adjacency Matrix
}

\description{
  Given a square binary adjacency matrix, the function
  \code{nbOrder} determines the integer matrix of neighbourhood orders
  (shortest-path distance).
}

\usage{
nbOrder(neighbourhood, maxlag = Inf)
}

\arguments{
  \item{neighbourhood}{
    a square, numeric or logical, and usually symmetric matrix with
    finite entries (and usually zeros on the diagonal) which indicates
    vertex adjacencies, i.e., first-order neighbourhood (interpreted as
    \code{neighbourhood == 1}, \emph{not} \code{>0}).
  }
  \item{maxlag}{
    positive scalar integer specifying an upper bound for the
    neighbourhood order. The default (\code{Inf}) means no truncation
    (but orders cannot be larger than the number of regions minus 1),
    whereas \code{maxlag = 1} just returns the input
    neighbourhood matrix (converted to binary integer mode).
  }
}

\value{
An integer matrix of neighbourhood orders, i.e., the shortest-path
distance matrix of the vertices.
The \code{dimnames} of the input \code{neighbourhood} matrix are preserved.
}

\author{
Sebastian Meyer
}

\seealso{
  \code{\link[spdep]{nblag}} from the \pkg{spdep} package
}

\examples{
## generate adjacency matrix
set.seed(1)
n <- 6
adjmat <- matrix(0, n, n)
adjmat[lower.tri(adjmat)] <- sample(0:1, n*(n-1)/2, replace=TRUE)
adjmat <- adjmat + t(adjmat)
adjmat

## determine neighbourhood order matrix
nblags <- nbOrder(adjmat)
nblags
}

\keyword{spatial}
\keyword{utilities}
