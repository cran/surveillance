\name{aggregate-methods}
\docType{methods}
\alias{aggregate.sts}
\alias{aggregate,sts-method}

\title{Aggregate an \code{"sts"} Object Over Time or Across Units}
\description{
  Aggregate the matrix slots of an \code{"\linkS4class{sts}"} object.
  Either the time series is aggregated so a new sampling frequency of
  \code{nfreq} observations per year is obtained (i.e., as in
  \code{\link{aggregate.ts}}), or the aggregation is over all
  columns (units).
}

\usage{
\S4method{aggregate}{sts}(x, by = "time", nfreq = "all", ...)
}

\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{by}{a string being either \code{"time"} or \code{"unit"}.}
  \item{nfreq}{new sampling frequency for \code{by="time"}. If
    \code{nfreq="all"} then all time points are summed.}
  \item{\dots}{unused (argument of the generic).}
}

\value{
  an object of class \code{"sts"}.
}

\section{Warning}{
  Aggregation over units fills the upperbound slot with
  \code{NA}s and the \code{map} slot is left as-is, but the object 
  cannot be plotted by unit any longer.
  
  The \code{populationFrac} slot is aggregated just like \code{observed}.
  Population fractions are recomputed if and only if \code{x} is no
  \code{multinomialTS} and already contains population fractions.
  This might not be intended, especially for aggregation over time.
}

\examples{
data("ha.sts")
dim(ha.sts)
dim(aggregate(ha.sts, by = "unit"))
dim(aggregate(ha.sts, nfreq = 13))
\dontshow{
## population(ha.sts) are trivial fractions, aggregate() should keep them
stopifnot(population(aggregate(ha.sts)) == 1/ncol(ha.sts))
## failed in surveillance <= 1.16.2
}
}

\keyword{methods}
