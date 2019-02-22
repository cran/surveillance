\name{epidata_intersperse}
\alias{intersperse}

\title{
Impute Blocks for Extra Stops in \code{"epidata"} Objects
}

\description{
  This function modifies an object inheriting from class \code{"epidata"} such 
  that it features the specified stop time points.  For this purpose, the time 
  interval in the event history into which the new stop falls will be split 
  up into two parts, one block for the time period until the new stop -- where 
  no infection or removal occurs -- and the other block for the time period 
  from the new stop to the end of the original interval.\cr
  Main application is to enable the use of \code{knots} in \code{twinSIR}, which
  are not existing stop time points in the \code{"epidata"} object.
}

\usage{
intersperse(epidata, stoptimes, verbose = FALSE)
}

\arguments{
  \item{epidata}{
    an object inheriting from class \code{"epidata"}.
  }
  \item{stoptimes}{
    a numeric vector of time points inside the observation period of the
    \code{epidata}.
  }
  \item{verbose}{
    logical indicating if a \code{\link{txtProgressBar}} should be shown
    while inserting blocks for extra \code{stoptimes}.
  }
}

\value{
  an object of the same class as \code{epidata} with additional time blocks
  for any new \code{stoptimes}.
}

\author{
Sebastian Meyer
}

\seealso{
  \code{\link{as.epidata.epidataCS}} where this function is used.
}

\examples{
data("hagelloch")
subset(hagelloch, start < 25 & stop > 25 & id \%in\% 9:13, select = 1:7)
# there is no "stop" time at 25, but we can add this extra stop
nrow(hagelloch)
moreStopsEpi <- intersperse(hagelloch, stoptimes = 25)
nrow(moreStopsEpi)
subset(moreStopsEpi, (stop == 25 | start == 25) & id \%in\% 9:13, select = 1:7)
}

\keyword{spatial}
\keyword{manip}
