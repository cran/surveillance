\name{[,sts-method}
\docType{methods}
\title{Subsetting \code{"sts"} Objects}
\alias{[,sts-method}
\alias{[,sts,ANY,ANY,ANY-method}
\description{
  The \code{"["}-method extracts parts of an
  \code{"\linkS4class{sts}"} object
  using row (time) and column (unit) indices.
}
\details{
  Column indices work as usual when indexing matrices,
  so may select units by name, position or a vector of booleans.
  Row indices can only be integer or logical. The \code{start} and
  \code{epoch} slots of the time series are adjusted accordingly.

  When subsetting columns, population fractions are recomputed if and
  only if \code{x} is no \code{multinomialTS} and already contains
  population fractions.
  
  A warning is issued if an irregular integer sequence is used to
  extract rows, e.g., \code{x[c(1,2,4),]}, which could destroy the
  structure of the time series (\code{freq}).
  
  \code{NA} indices are not supported, negative indices are.

  The formal \code{drop} argument is ignored (always \code{FALSE}).

  Note that a \code{[<-} method (i.e., subassignment) is not implemented.
}

\examples{
data("ha.sts")
haagg <- aggregate(ha.sts, nfreq=13)

plot(haagg[, 3])       # Single series
plot(haagg[1:30, 3])   # Somewhat shorter

#Counts at time 20
plot(haagg[20, ], type = observed ~ unit)
}

\keyword{methods}
\keyword{array}
