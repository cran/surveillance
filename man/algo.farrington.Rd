\name{algo.farrington}
\alias{algo.farrington}
\alias{farrington}
\encoding{latin1}

\title{Surveillance for Count Time Series Using the Classic Farrington Method}
\description{
  Implements the procedure of Farrington et al. (1996).
  At each time point of the specified \code{range}, a GLM is fitted to
  predict the counts. This is then compared to the observed
  counts. If the observation is above a specific quantile of
  the prediction interval, then an alarm is raised.
}
\usage{
# original interface for a single "disProg" time series
algo.farrington(disProgObj, control=list(
    range=NULL, b=5, w=3, reweight=TRUE, verbose=FALSE, plot=FALSE,
    alpha=0.05, trend=TRUE, limit54=c(5,4), powertrans="2/3",
    fitFun="algo.farrington.fitGLM.fast"))

# wrapper for "sts" data, possibly multivariate
farrington(sts, control=list(
    range=NULL, b=5, w=3, reweight=TRUE, verbose=FALSE,
    alpha=0.05), ...)
}
\arguments{
  \item{disProgObj}{
    an object of class \code{"disProg"} (a list including \code{observed} and
    \code{state} time series).
  }
  \item{control}{list of control parameters
    \describe{
    \item{\code{range}}{Specifies the index of all timepoints which
        should be tested. If range is \code{NULL} the maximum number
        of possible weeks is used (i.e. as many weeks as possible while
	still having enough reference values).}
    \item{\code{b}}{how many years back in time to include when
        forming the base counts.}
    \item{\code{w}}{windows size, i.e. number of weeks to include
        before and after the current week}
    \item{\code{reweight}}{Boolean specifying whether to perform reweight step}
    \item{\code{trend}}{If \code{TRUE} a trend is included and kept in
        case the conditions documented in Farrington et al. (1996) are met
        (see the results). If \code{FALSE} then NO trend is fit.}
    \item{\code{verbose}}{Boolean indicating whether to show extra debugging information.}
    \item{\code{plot}}{Boolean specifying whether to show the final GLM model fit graphically (use
        History|Recording to see all pictures).}
    \item{\code{powertrans}}{Power transformation to apply to the
        data. Use either "2/3" for skewness correction (Default),
        "1/2" for variance stabilizing transformation or "none" for no
        transformation.}
    \item{\code{alpha}}{An approximate (two-sided) \eqn{(1-\alpha)}
        prediction interval is calculated.}
    \item{\code{limit54}}{To avoid alarms in cases where the time series only
        has about 0-2 cases the algorithm uses the following heuristic
        criterion (see Section 3.8 of the Farrington paper) to protect
        against low counts: no alarm is sounded if fewer than
        \eqn{cases=5} reports were received in the past \eqn{period=4}
        weeks. \code{limit54=c(cases,period)} is a vector allowing the
        user to change these numbers. Note: As of version 0.9-7 the
        term "last" period of weeks includes the current week -
        otherwise no alarm is sounded for horrible large numbers if
        the four weeks before that are too low.}
    \item{\code{fitFun}}{String containing the name of the fit
        function to be used for fitting the GLM. The options are
        \code{algo.farrington.fitGLM.fast} (default) and
        \code{algo.farrington.fitGLM} or 
        \code{algo.farrington.fitGLM.populationOffset}. See details of
        \code{\link{algo.farrington.fitGLM}} for more information.}
    }
  }
  \item{sts}{an object of class \code{"\link{sts}"}.}
  \item{\dots}{arguments for \code{\link{wrap.algo}},
    e.g., \code{verbose=FALSE}.}
}
\details{
  The following steps are performed according to the Farrington
  et al. (1996) paper.
\enumerate{
\item fit of the initial model and initial estimation of mean and
  overdispersion.
\item calculation of the weights omega (correction for past outbreaks)
\item refitting of the model
\item revised estimation of overdispersion
\item rescaled model
\item omission of the trend, if it is not significant
\item repetition of the whole procedure
\item calculation of the threshold value
\item computation of exceedance score
}
}
\value{
For \code{algo.farrington}, a list object of class \code{"survRes"}
with elements \code{alarm}, \code{upperbound}, \code{trend},
\code{disProgObj}, and \code{control}.

For \code{farrington}, the input \code{"\link{sts}"} object with updated
\code{alarm}, \code{upperbound} and \code{control} slots, and subsetted
to \code{control$range}.
}

\examples{
#load "disProg" data
data("salmonella.agona")

#Do surveillance for the last 42 weeks
n <- length(salmonella.agona$observed)
control <- list(b=4,w=3,range=(n-42):n,reweight=TRUE, verbose=FALSE,alpha=0.01)
res <- algo.farrington(salmonella.agona,control=control)
plot(res)

#Generate Poisson counts and create an "sts" object
set.seed(123)
x <- rpois(520,lambda=1)
stsObj <- sts(observed=x, frequency=52)

if (surveillance.options("allExamples")) {
#Compare timing of the two possible fitters for algo.farrington
  range <- 312:520
  system.time( sts1 <- farrington(stsObj, control=list(range=range,
                         fitFun="algo.farrington.fitGLM.fast"), verbose=FALSE))
  system.time( sts2 <- farrington(stsObj, control=list(range=range,
                         fitFun="algo.farrington.fitGLM"), verbose=FALSE))
  #Check if results are the same
  stopifnot(upperbound(sts1) == upperbound(sts2))
}
}

\author{M. \enc{Höhle}{Hoehle}}
\seealso{
  \code{\link{algo.farrington.fitGLM}},
  \code{\link{algo.farrington.threshold}}

  An improved Farrington algorithm is available as function
  \code{\link{farringtonFlexible}}.
}
\keyword{classif}
\references{
  A statistical algorithm for the early detection of outbreaks of
  infectious disease, Farrington, C.P., Andrews, N.J, Beale A.D. and
  Catchpole, M.A. (1996), J. R. Statist. Soc. A, 159, 547-563.
}
