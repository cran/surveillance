\name{refvalIdxByDate}
\alias{refvalIdxByDate}
\title{Compute indices of reference value using Date class}
\description{ 
  The reference values are formed based on computations 
  of \code{seq} for Date class arguments.  
}
\usage{
refvalIdxByDate(t0, b, w, epochStr, epochs)
}
\arguments{
\item{t0}{A Date object describing the time point}
\item{b}{Number of years to go back in time}
\item{w}{Half width of window to include reference values for}
\item{epochStr}{One of \code{"1 month"}, \code{"1 week"} or \code{"1 day"}}
\item{epochs}{Vector containing the epoch value of the sts/disProg object}
}
\details{ 
  Using the Date class the reference values are formed as follows:
  Starting from \code{t0} go i, i= 1,...,\code{b} years back in time.
  For each year, go \code{w} epochs back and include from here to
  \code{w} epochs after \code{t0}.

  In case of weeks we always go back to the closest Monday of this
  date. In case of months we also go back in time to closest 1st of
  month. 
}
\value{
  a vector of indices in epochs which match
}
\keyword{chron}
