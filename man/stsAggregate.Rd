\name{aggregate-methods}
\docType{methods}
\alias{aggregate,sts,ANY,ANY-method}
\alias{aggregate,sts-method}

\title{Aggregate an \code{"sts"} Object Over Time or Across Units}
\description{
  Method to aggregate the matrix slots of an sts object.
  Either the time series is aggregated so a new sampling frequency of
  \code{nfreq} units per time slot is obtained (i.e as in
  \code{\link{aggregate.ts}}) or the aggregation is over all
  \code{ncol} units.

  Note: The function is not 100\% consistent with what the generic
  function \code{\link{aggregate}} does.
}

\section{Methods}{
  \describe{
    \item{x = "sts", by="time", nfreq="all",...}{
      \describe{
	\item{x}{an object of class \code{sts}}
	\item{by}{a string being either "time" or "unit"}
	\item{nfreq}{new sampling frequency if \code{by=="time"}. If
	  \code{nfreq=="all"} then all time instances are summed.}
	\item{...}{not used}
      }
      returns an object of class \code{sts}
    }
  }
}

\section{Warning}{
  Aggregation by unit sets the upperbound slot to
  \code{NA} and the \code{map} slot is left as-is, but the object 
  cannot be plotted by unit any longer. 
}
  
\seealso{
  \code{\link{aggregate}}
}

\examples{
  data("ha.sts")
  dim(ha.sts)
  dim(aggregate(ha.sts,by="unit"))
  dim(aggregate(ha.sts,nfreq=13))
}

\keyword{methods}
