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
\S4method{[}{sts}(x, i, j, ..., drop = FALSE)
}
\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{i}{optional row index (integer or logical vector).}
  \item{j}{optional column index (character, integer, or logical vector).}
  \item{drop}{logical: Should subsetting by \code{j} be applied to
    the map as well? This requires a \emph{character} index and is
    disabled by default. It does not affect the other slots:
    Dimensions are never dropped.}
  \item{\dots}{ignored.}
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

# Show a (subset of a) single time series
plot(ha.sts[,7])
plot(ha.sts[year(ha.sts)==2006, 7])

# Map a single time point
plot(ha.sts[5*52+26,], type=observed~unit)
plot(ha.sts,           type=observed~unit, tps=5*52+26) # same -> ?stsplot_space

# Restrict the data (and the map) to a subset of the districts
plot(ha.sts[,c("pank","lich")],           type=observed~unit, labels=TRUE)
plot(ha.sts[,c("pank","lich"),drop=TRUE], type=observed~unit, labels=TRUE)
}

\keyword{methods}
