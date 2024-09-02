\name{algo.rki}
\alias{algo.rkiLatestTimepoint}
\alias{algo.rki}
\alias{algo.rki1}
\alias{algo.rki2}
\alias{algo.rki3}
\encoding{latin1}
\title{The system used at the RKI}
\description{
  Evaluation of timepoints with the detection algorithms used by the RKI
}
\usage{
algo.rkiLatestTimepoint(disProgObj, timePoint = NULL,
                        control = list(b = 2, w = 4, actY = FALSE))
algo.rki(disProgObj, control = list(range = range,
         b = 2, w = 4, actY = FALSE))
algo.rki1(disProgObj, control = list(range = range))
algo.rki2(disProgObj, control = list(range = range))
algo.rki3(disProgObj, control = list(range = range))
}
\arguments{
\item{disProgObj}{object of class disProg (including the observed and the state chain).}
\item{timePoint}{time point which should be evaluated in \code{algo.rkiLatestTimepoint}. The
    default is to use the latest timepoint.}
\item{control}{control object: \code{range} determines the desired
    timepoints which should be evaluated, \code{b} describes the number of years to go
    back for the reference values, \code{w} is the half window width for the reference
    values around the appropriate timepoint and \code{actY} is a boolean to decide if
    the year of \code{timePoint} also spend \code{w} reference values of the past.
    As default \code{b}, \code{w}, \code{actY} are set for the RKI 3 system.
  }
}
\value{
  \code{algo.rkiLatestTimepoint} returns a list of class \code{survRes} (surveillance result), which
  includes the alarm value (alarm = 1, no alarm = 0) for recognizing an
  outbreak, the threshold value for recognizing the alarm and
  the input object of class disProg.
  
  \code{algo.rki} gives a list of class \code{survRes} which includes the vector
  of alarm values for every timepoint in \code{range}, the vector of threshold values
  for every timepoint in \code{range} for the system specified by \code{b}, \code{w} and
  \code{actY}, the range and the input object of class disProg.
  \code{algo.rki1} returns the same for the RKI 1 system, \code{algo.rki2}
  for the RKI 2 system and \code{algo.rki3} for the RKI 3 system.
}
\details{
  Using the reference values for calculating an upper limit (threshold),
  alarm is given if the actual value is bigger than a computed threshold.
  \code{algo.rki} calls \code{algo.rkiLatestTimepoint} for the values specified
  in \code{range} and for the system specified in \code{control}.
  \code{algo.rki1} calls \code{algo.rkiLatestTimepoint} for the values specified
  in \code{range} for the RKI 1 system.
  \code{algo.rki2} calls \code{algo.rkiLatestTimepoint} for the values specified
  in \code{range} for the RKI 2 system.
  \code{algo.rki3} calls \code{algo.rkiLatestTimepoint} for the values specified
  in \code{range} for the RKI 3 system.
  \itemize{
  \item \code{"RKI 1"} reference values from 6 weeks ago
  \item \code{"RKI 2"} reference values from 6 weeks ago and
    13 weeks of the year ago (symmetrical around the
    comparable week).
  \item \code{"RKI 3"}  18 reference values. 9 from the year ago
    and 9 from two years ago (also symmetrical around the
    comparable week).
  }
}
\seealso{
  \code{\link{algo.bayesLatestTimepoint}} and \code{\link{algo.bayes}} for
  the Bayes system.
}
\author{M. \enc{Höhle}{Hoehle}, A. Riebler, Christian Lang}
\examples{
# Create a test object
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208, A = 1,
                              alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.7)

# Test week 200 to 208 for outbreaks with a selfdefined rki
algo.rki(disProgObj, control = list(range = 200:208, b = 1,
                                    w = 5, actY = TRUE))
# The same for rki 1 to rki 3
algo.rki1(disProgObj, control = list(range = 200:208))
algo.rki2(disProgObj, control = list(range = 200:208))
algo.rki3(disProgObj, control = list(range = 200:208))

# Test for rki 1 the latest timepoint
algo.rkiLatestTimepoint(disProgObj)
}
\keyword{classif}
