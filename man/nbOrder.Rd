\name{nbOrder}
\alias{nbOrder}
\title{
  Determine Neighbourhood Order Matrix from Binary Adjacency Matrix
}

\description{
  Given a square binary adjacency matrix, the function
  \code{nbOrder} determines the integer matrix of neighbourhood orders
  (shortest-path distance) using the function \code{\link[spdep]{nblag}}
  from the \pkg{spdep} package.
}

\usage{
nbOrder(neighbourhood, maxlag = 1)
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
    neighbourhood order. The default (1) just returns the input
    neighbourhood matrix (converted to binary integer mode).
    \code{maxlag} is automatically trimmed to one less than the number
    of regions (there cannot be higher orders) and then converted to
    integer, thus, \code{maxlag = Inf} also works.
  }
}

\value{
An integer matrix of neighbourhood orders, i.e., the shortest-path
distance matrix of the vertices.
The \code{dimnames} of the input \code{neighbourhood} matrix are preserved.
}

\note{
By the end, the function issues a \code{\link{message}} informing about the
range of maximum neighbourhood order by region.
}

\author{
Sebastian Meyer
}

\seealso{
  \code{\link[spdep]{nblag}} from the \pkg{spdep} package, on which this
  wrapper depends.
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
if (requireNamespace("spdep")) {
    nbmat <- nbOrder(adjmat, maxlag=Inf)
    nbmat
}
}

\keyword{spatial}
\keyword{utilities}
