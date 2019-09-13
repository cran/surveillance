\name{stsXtrct}
\docType{methods}
\title{Subsetting \code{"sts"} Objects}
\alias{[,sts-method} % for convenience
\alias{[,sts,ANY,ANY,ANY-method}
\description{
  The \code{[}-method extracts parts of an
  \code{"\linkS4class{sts}"} object
  using row (time) and column (unit) indices.
}
\usage{
\S4method{[}{sts}(x, i, j, ..., drop)
}
\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{i}{row index (integer or logical vector).}
  \item{j}{column index (character, integer, or logical vector).}
  \item{\dots,drop}{unused (arguments of the generic).\cr
    Dimensions are never dropped.}
}
\value{
  an object of class \code{"sts"}.
}
\details{
  Row indices are used to select a subset of the original time period.
  The \code{start} and \code{epoch} slots of the time series are
  adjusted accordingly.
  A warning is issued if an irregular integer sequence is used to
  extract rows, e.g., \code{x[c(1,2,4),]}, which could destroy the
  structure of the time series (\code{freq}).

  Column indices work as usual when indexing matrices,
  so may select units by name, position or a vector of booleans.
  When subsetting columns, population fractions are recomputed if and
  only if \code{x} is no \code{multinomialTS} and already contains
  population fractions.
    
  \code{NA} indices are not supported, negative indices are.

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
