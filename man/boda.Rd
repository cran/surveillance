\encoding{latin1}
\name{boda}
\alias{boda}

\title{Bayesian Outbreak Detection Algorithm (BODA)}

\description{
  The function takes \code{range} values of a univariate surveillance time
  series \code{sts} and for each time point uses a negative binomial
  regression model to compute the predictive posterior distribution for
  the current observation. The
  \eqn{(1-\alpha)\cdot 100\%}{(1-alpha)*100\%}
  quantile of this predictive distribution is then
  used as bound: If the actual observation is above the bound an alarm
  is raised. 
  The Bayesian Outbreak Detection Algorithm (\code{boda}) is due to
  Manitz and \enc{Höhle}{Hoehle} (2013) and its implementation is
  illustrated in Salmon et al. (2016).
  However, \code{boda} should be considered as an experiment, see the
  Warning section below!
}

\usage{
boda(sts, control = list(
    range=NULL, X=NULL, trend=FALSE, season=FALSE,
    prior=c('iid','rw1','rw2'), alpha=0.05, mc.munu=100, 
    mc.y=10, verbose=FALSE,multicore=TRUE,
    samplingMethod=c('joint','marginals'),
    quantileMethod=c("MC","MM")
))
}

\arguments{
  \item{sts}{object of class sts (including the \code{observed} and the
    \code{state} time series)}
  \item{control}{Control object given as a \code{list} containing the following components:
    \describe{
      \item{\code{range}}{Specifies the index of all timepoints which
	should be tested. If range is \code{NULL} all possible
	timepoints are used.} 
      \item{\code{X}}{}
      \item{\code{trend}}{Boolean indicating whether a linear trend term should be
	included in the model for the expectation the log-scale}
      \item{\code{season}}{Boolean to indicate whether a cyclic spline
	should be included.}
      \item{\code{alpha}}{The threshold for declaring an observed count as
	an aberration is the \eqn{(1-\alpha)\cdot 100\%}{(1-alpha)*100\%}
	quantile of the predictive posterior.}
      \item{\code{mc.munu}}{}
      \item{\code{mc.y}}{Number of samples of \eqn{y}{y} to generate for
	each par of the mean and size parameter. A total of \eqn{mc.munu
	  \times mc.y}{mc.munu*mc.y} samples are generated.}
      \item{\code{verbose}}{Argument sent to the inla call. When using ESS
	it might be necessary to force verbose mode for INLA to work.}
      \item{\code{multicore}}{Detect using \code{parallel::detectCores}
	how many logical cores are available and set INLA to use this number.}
      \item{\code{samplingMethod}}{Should one sample from the parameters joint distribution (joint) or from their respective marginal posterior distribution (marginals)?}
      \item{quantileMethod}{Character, either \code{MC} or \code{MM}. Indicates how to compute the quantile based on the posterior distribution (no matter the inference method): either by sampling \code{mc.munu} values from the posterior distribution of the parameters and then for each sampled parameters vector sampling \code{mc.y} response values so that one gets a vector of response values based on which one computes an empirical quantile (MC method, as explained in Manitz and \enc{Höhle}{Hoehle} 2013); or by sampling \code{mc.munu} from the posterior distribution of the parameters and then compute the quantile of the mixture distribution using bisectioning, which is faster.}
    }
  }
}

\note{
  This function requires the \R package \pkg{INLA}, which is currently
  \emph{not} available from CRAN. It can be obtained from INLA's own
  repository via
  \code{install.packages("INLA", repos="https://inla.r-inla-download.org/R/stable")}.
}

\section{Warning}{
  This function is currently experimental!! It also heavily
  depends on the \pkg{INLA} package so changes there might affect the
  operational ability of this function. Since the computations
  for the Bayesian GAM are quite involved do not expect this function
  to be particularly fast. Future work could focus on improving the
  speed, e.g., one issue would be to make the inference work in a sequential
  fashion.
}

\keyword{classif}

\examples{
\dontrun{
  ## running this example takes a couple of minutes

  #Load the campylobacteriosis data for Germany
  data("campyDE")
  #Make an sts object from the data.frame
  cam.sts <-  sts(epoch=campyDE$date,
                  observed=campyDE$case, state=campyDE$state)

  #Define monitoring period
#  range <- which(epoch(cam.sts)>=as.Date("2007-01-01"))
#  range <- which(epoch(cam.sts)>=as.Date("2011-12-10"))
  range <- tail(1:nrow(cam.sts),n=2)

  control <- list(range=range, X=NULL, trend=TRUE, season=TRUE,
                  prior='iid', alpha=0.025, mc.munu=100, mc.y=10,
                  samplingMethod = "joint")

  #Apply the boda algorithm in its simples form, i.e. spline is
  #described by iid random effects and no extra covariates
  library("INLA")  # needs to be attached
  cam.boda1 <- boda(cam.sts, control=control)

  plot(cam.boda1, xlab='time [weeks]', ylab='No. reported', dx.upperbound=0)
}
}

\author{J. Manitz, M. \enc{Höhle}{Hoehle}, M. Salmon}

\references{
  Manitz, J. and \enc{Höhle}{Hoehle}, M. (2013):
  Bayesian outbreak detection algorithm for monitoring reported cases of
  campylobacteriosis in Germany.
  Biometrical Journal, 55(4), 509-526.

  Salmon, M., Schumacher, D. and \enc{Höhle}{Hoehle}, M. (2016):
  Monitoring count time series in \R: Aberration detection in public
  health surveillance. \emph{Journal of Statistical Software},
  \bold{70} (10), 1-35. \doi{10.18637/jss.v070.i10}
}
