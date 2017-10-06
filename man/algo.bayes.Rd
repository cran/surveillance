\name{algo.bayes}
\alias{algo.bayes}
\alias{algo.bayesLatestTimepoint}
\alias{algo.bayes1}
\alias{algo.bayes2}
\alias{algo.bayes3}
\encoding{latin1}

\title{The Bayes System}
\description{
        Evaluation of timepoints with the Bayes subsystem
        1, 2, 3 or a self defined Bayes subsystem.
}
\usage{
  algo.bayesLatestTimepoint(disProgObj, timePoint = NULL,
       control = list(b = 0, w = 6, actY = TRUE,alpha=0.05))
  algo.bayes(disProgObj, control = list(range = range,
       b = 0, w = 6, actY = TRUE,alpha=0.05))
  algo.bayes1(disProgObj, control = list(range = range))
  algo.bayes2(disProgObj, control = list(range = range))
  algo.bayes3(disProgObj, control = list(range = range))
}
\arguments{
        \item{disProgObj}{object of class disProg (including the observed and the state chain)}
        \item{timePoint}{time point which should be evaluated in
	  \code{algo.bayes LatestTimepoint}. The
                default is to use the latest timepoint}
        \item{control}{control object: \code{range} determines the desired
                timepoints which should be evaluated, \code{b} describes the number of years to go
                back for the reference values, \code{w} is the half window width for the reference
                values around the appropriate timepoint and \code{actY} is a boolean to decide if
                the year of \code{timePoint} also contributes \code{w} reference values. The parameter \code{alpha} is the \eqn{(1-\alpha)}-quantile to use in order to calculate the upper threshold.
                As default \code{b}, \code{w}, \code{actY} are set for the Bayes 1 system with \code{alpha}=0.05.
        }
}
\value{
\item{survRes}{
    \code{algo.bayesLatestTimepoint} returns a list of class \code{survRes} (surveillance result), which
    includes the alarm value for recognizing an
    outbreak (1 for alarm, 0 for no alarm), the threshold value for recognizing the alarm and
    the input object of class disProg.
    \code{algo.bayes} gives a list of class \code{survRes} which includes the vector
    of alarm values for every timepoint in \code{range} and the vector of threshold values
    for every timepoint in \code{range} for the system specified by \code{b}, \code{w} and
    \code{actY}, the range and the input object of class disProg.
    \code{algo.bayes1} returns the same for the Bayes 1 system, \code{algo.bayes2}
    for the Bayes 2 system and \code{algo.bayes3} for the Bayes 3 system.
  }
}
\details{
  Using the reference values the \eqn{(1-\alpha)\cdot
    100\%}{(1-alpha)*100\%} quantile of the
  predictive posterior distribution is calculated as a threshold.
  An alarm is given if the actual value is bigger or equal than this threshold.
  It is possible to show using analytical computations that the predictive
  posterior in this case is the negative
  binomial distribution. Note: \code{algo.rki} or \code{algo.farrington}
  use two-sided prediction intervals -- if one wants to compare with
  these procedures it is necessary to use an alpha, which is half the
  one used for these procedures.

  Note also that \code{algo.bayes} calls
  \code{algo.bayesLatestTimepoint} for the values specified in
  \code{range} and for the system specified in \code{control}.
  \code{algo.bayes1}, \code{algo.bayes2}, \code{algo.bayes3} call
  \code{algo.bayesLatestTimepoint} for the values specified in
  \code{range} for the Bayes 1 system, Bayes 2 system or Bayes 3 system.

  \itemize{
  \item \code{"Bayes 1"} reference values from 6 weeks. Alpha is fixed a
 t 0.05.
  \item \code{"Bayes 2"} reference values from 6 weeks ago and
    13 weeks of the previous year (symmetrical around the
    same week as the current one in the previous year). Alpha is fixed at 0.05.
  \item \code{"Bayes 3"}  18 reference values. 9 from the year ago
    and 9 from two years ago (also symmetrical around the
    comparable week). Alpha is fixed at 0.05.
  }

  The procedure is now able to handle \code{NA}'s in the reference
  values. In the summation and when counting the number of observed
  reference values these are simply not counted.
}
\seealso{
        \code{\link{algo.call}}, \code{\link{algo.rkiLatestTimepoint}} and \code{\link{algo.rki}} for
        the RKI system.
}
\author{M. \enc{Höhle}{Hoehle}, A. Riebler, C. Lang}
\examples{
    disProg <- sim.pointSource(p = 0.99, r = 0.5, length = 208, A = 1,
                                    alpha = 1, beta = 0, phi = 0,
                                    frequency = 1, state = NULL, K = 1.7)

    # Test for bayes 1 the latest timepoint
    algo.bayesLatestTimepoint(disProg)

    # Test week 200 to 208 for outbreaks with a selfdefined bayes
    algo.bayes(disProg, control = list(range = 200:208, b = 1,
                                                w = 5, actY = TRUE,alpha=0.05))
    # The same for bayes 1 to bayes 3
    algo.bayes1(disProg, control = list(range = 200:208,alpha=0.05))
    algo.bayes2(disProg, control = list(range = 200:208,alpha=0.05))
    algo.bayes3(disProg, control = list(range = 200:208,alpha=0.05))
}

\keyword{classif}
\source{
Riebler, A. (2004), Empirischer Vergleich von statistischen Methoden zur
Ausbruchserkennung bei Surveillance Daten, Bachelor's thesis.
}
