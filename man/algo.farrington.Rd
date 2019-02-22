\name{algo.farrington}
\alias{algo.farrington}
\encoding{latin1}

\title{Surveillance for Count Time Series Using the Classic Farrington Method}
\description{
  The function takes \code{range} values of the surveillance time
  series \code{disProgObj} and for each time point uses a GLM to
  predict the number of counts according to the procedure by
  Farrington et al. (1996). This is then compared to the observed
  number of counts. If the observation is above a specific quantile of
  the prediction interval, then an alarm is raised.
}
\usage{
algo.farrington(disProgObj, control=list(
    range=NULL, b=5, w=3, reweight=TRUE, verbose=FALSE, plot=FALSE,
    alpha=0.05, trend=TRUE, limit54=c(5,4), powertrans="2/3",
    fitFun="algo.farrington.fitGLM.fast"))
}
\arguments{
  \item{disProgObj}{
    object of class disProgObj (including the \code{observed} and the
    \code{state} time series.)
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
An object of class \code{"survRes"}.
}

\examples{
#Read Salmonella Agona data
data("salmonella.agona")

#Do surveillance for the last 100 weeks.
n <- length(salmonella.agona$observed)
#Set control parameters.
control <- list(b=4,w=3,range=(n-100):n,reweight=TRUE, verbose=FALSE,alpha=0.01)
res <- algo.farrington(salmonella.agona,control=control)
#Plot the result.
plot(res,disease="Salmonella Agona",method="Farrington")

\dontrun{
#Generate Poisson counts and convert into an "sts" object
set.seed(123)
x <- rpois(520,lambda=1)
sts <- sts(observed=x, state=x*0, freq=52)

#Compare timing of the two possible fitters for algo.farrington (here using S4)
system.time( sts1 <- farrington(sts, control=list(range=312:520,
                       fitFun="algo.farrington.fitGLM.fast")))
system.time( sts2 <- farrington(sts, control=list(range=312:520,
                       fitFun="algo.farrington.fitGLM")))

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
