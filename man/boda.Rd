\name{boda}
\alias{boda}
\encoding{latin1}

\title{Surveillance for an univariate count data time series using the
  Bayesian Outbreak Detection Algorithm (BODA) described in Manitz and
  \enc{Höhle}{Hoehle} (2013)}
\description{
%
  The function takes \code{range} values of the surveillance time
  series \code{sts} and for each time point uses a negative binomial
  regression model to compute the predictive posterior distribution for
  the current observation. The \eqn{(1-\alpha)\cdot
  100\%}{(1-alpha)*100} quantile of this predictive distribution is then
  used as bound: If the actual observation is above the bound an alarm
  is raised. 
%
}
\usage{
boda(sts, control=list(range=NULL, X=NULL, trend=FALSE,
     season=FALSE, prior=c('iid','rw1','rw2'), alpha=0.05,
     mc.munu=100, mc.y=10, verbose=FALSE,multicore=TRUE))
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
      \item{\code{season}}{Boolean to indicate wheather a cyclic spline
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
    }
  }
}
\details{

  Note: This function requires presence of the INLA R package, which is
  NOT available from CRAN. It can can be downloaded by calling
  \code{source("http://www.math.ntnu.no/inla/givemeINLA.R")} as described in
  detail at \url{http://www.r-inla.org/download}. 

  WARNING: This function is currently experimental!! It also heavily
  depends on the INLA package so changes here might affect the
  operational ability of the function. Since the computations
  for the Bayesian GAM are quite involved do not expected this function
  to be particularly fast. Future work could focus on improving the
  speed, e.g. one issue would be to make the inference work in sequential
  fashion.
}

\keyword{classif}

\examples{
  #Load the campylobacteriosis data for Germany
  data("campyDE")
  #Make an sts object from the data.frame
  cam.sts <-  new("sts",epoch=as.numeric(campyDE$date),observed=campyDE$case,
                   state=campyDE$state, epochAsDate=TRUE)

  \dontrun{
  #Trial run code. To be made an actual working example.
  #source("../R/boda.R")

  # define monitoring period
#  range <- which(epoch(cam.sts)>=as.Date("2007-01-01"))
#  range <- which(epoch(cam.sts)>=as.Date("2011-12-10"))
  range <- tail(1:nrow(cam.sts),n=2)

  control <- list(range=range, X=NULL, trend=TRUE, season=TRUE,
                  prior='iid', alpha=0.025, mc.munu=100, mc.y=10)

  #Apply the boda algorithm in its simples form, i.e. spline is
  #described by iid random effects and no extra covariates
  #(NOTE: requires the INLA package to be installed)
  cam.boda1 <- boda(cam.sts, control=control)

  #In case INLA is not installed, boda throws an error
  if(!inherits(cam.boda1,'try-error')){
    plot(cam.boda1,xlab='time [weeks]', ylab='No. reported',dx.upperbound=0)
  }

  } 
}

\author{J. Manitz and M. \enc{Höhle}{Hoehle}}

\references{
  Bayesian model algorithm for monitoring reported cases of
  campylobacteriosis in Germany (2013), Manitz J and Höhle M,
  Biometrical Journal, 55(4), pp. 509 526.
}


