\encoding{latin1}
\name{bodaDelay}
\alias{bodaDelay}
\title{Bayesian Outbreak Detection in the Presence of Reporting Delays}
\usage{
bodaDelay(sts, control = list(
  range = NULL, b = 5, w = 3, mc.munu = 100, mc.y = 10,
  pastAberrations = TRUE, verbose = FALSE,
  alpha = 0.05, trend = TRUE, limit54 = c(5,4), 
  inferenceMethod = c("asym","INLA"), quantileMethod = c("MC","MM"),
  noPeriods = 1, pastWeeksNotIncluded = NULL, delay = FALSE))
}

\arguments{
\item{sts}{sts-object to be analysed. Needs to have a reporting triangle.}

\item{control}{list of control arguments:
\describe{
\item{\code{b}}{How many years back in time to include when forming the base counts.}

\item{\code{w}}{Window's half-size, i.e. number of weeks to include before and after the current week in each year.}

\item{\code{range}}{Specifies the index of all timepoints which should be tested. If range is \code{NULL} all possible timepoints are used.}

\item{\code{pastAberrations}}{Boolean indicating whether to include an effect for past outbreaks
in a second fit of the model. This option only makes sense if \code{inferenceMethod} is \code{INLA}, as it is not supported by the other inference method.}

\item{\code{verbose}}{Boolean specifying whether to show extra debugging information.}

\item{\code{alpha}}{An approximate (one-sided) \eqn{(1-\alpha)\cdot 100\%}
prediction interval is calculated unlike the original method where it was a two-sided interval. The upper limit of this interval
i.e. the \eqn{(1-\alpha)\cdot 100\%} quantile serves as an upperbound.}

\item{\code{trend}}{Boolean indicating whether a trend should be included}

\item{\code{noPeriods}}{Number of levels in the factor allowing to use more baseline. If
equal to 1 no factor variable is created, the set of reference values is defined as in
Farrington et al (1996).}

\item{\code{inferenceMethod}}{Which inference method used, as defined in Salmon et al. (2015). If one chooses \code{"INLA"} then inference is performed with INLA. If one chooses \code{"asym"} (default) then the asymptotic normal approximation of the posteriori is used.}

\item{\code{pastWeeksNotIncluded}}{Number of past weeks to ignore in the calculation.
The default (\code{NULL}) means to use the value of \code{control$w}.}

\item{\code{delay}}{Boolean indicating whether to take reporting delays into account.}

\item{\code{mc.munu}}{Number of samples for the parameters of the negative binomial distribution for calculating a threshold}

\item{\code{mc.y}}{Number of samples for observations
when performing Monte Carlo to calculate a threshold}

\item{\code{limit54}}{c(cases,period) is a vector allowing the
user to change these numbers.}

\item{\code{quantileMethod}}{Character, either \code{"MC"} (default) or \code{"MM"}. Indicates how to compute the quantile based on the posterior distribution (no matter the inference method): either by sampling \code{mc.munu} values from the posterior distribution of the parameters and then for each sampled parameters vector sampling \code{mc.y} response values so that one gets a vector of response values based on which one computes an empirical quantile (MC method, as explained in Salmon et al. 2015); or by sampling \code{mc.munu} from the posterior distribution of the parameters and then compute the quantile of the mixture distribution using bisectioning, which is faster.}
}
}
}

\description{
The function takes \code{range} values of the surveillance time
series \code{sts} and for each time point uses a Bayesian model of the negative binomial family with
log link inspired by the work of Noufaily et al. (2012) and of Manitz and \enc{Höhle}{Hoehle} (2014). It allows delay-corrected aberration detection as explained in Salmon et al. (2015). A \code{reportingTriangle} has to be provided in the \code{control} slot.
}

\examples{
\dontrun{
data("stsNewport")
salm.Normal <- list()
salmDelayAsym <- list()
for (week in 43:45){
  listWeeks <- as.Date(row.names(stsNewport@control$reportingTriangle$n))
  dateObs <- listWeeks[isoWeekYear(listWeeks)$ISOYear==2011 &
                       isoWeekYear(listWeeks)$ISOWeek==week]
  stsC <- sts_observation(stsNewport,
                          dateObservation=dateObs,
                          cut=TRUE)
  inWeeks <- with(isoWeekYear(epoch(stsC)),
                  ISOYear == 2011 & ISOWeek >= 40 & ISOWeek <= 48)
  
  rangeTest <- which(inWeeks)
  alpha <- 0.07

  # Control slot for Noufaily method          
  controlNoufaily <- list(range=rangeTest,noPeriods=10,
                          b=4,w=3,weightsThreshold=2.58,pastWeeksNotIncluded=26,
                          pThresholdTrend=1,thresholdMethod="nbPlugin",alpha=alpha*2,
                          limit54=c(0,50))
  
  # Control slot for the Proposed algorithm with D=0 correction
  controlNormal <- list(range = rangeTest, b = 4, w = 3,
                        reweight = TRUE, mc.munu=10000, mc.y=100,
                        verbose = FALSE,
                        alpha = alpha, trend = TRUE,
                        limit54=c(0,50), 
                        noPeriods = 10, pastWeeksNotIncluded = 26,
                        delay=FALSE)
  
  # Control slot for the Proposed algorithm with D=10 correction
  controlDelayNorm <-  list(range = rangeTest, b = 4, w = 3,
                            reweight = FALSE, mc.munu=10000, mc.y=100,
                            verbose = FALSE,
                            alpha = alpha, trend = TRUE,
                            limit54=c(0,50), 
                            noPeriods = 10, pastWeeksNotIncluded = 26,
                            delay=TRUE,inferenceMethod="asym")
  
  set.seed(1)
  salm.Normal[[week]] <- farringtonFlexible(stsC, controlNoufaily)
  salmDelayAsym[[week]] <- bodaDelay(stsC, controlDelayNorm)
}

opar <- par(mfrow=c(2,3))
lapply(salmDelayAsym[c(43,44,45)],plot, legend=NULL, main="", ylim=c(0,35))
lapply(salm.Normal[c(43,44,45)],plot, legend=NULL, main="", ylim=c(0,35))
par(opar)
}
}

\references{
  Farrington, C.P., Andrews, N.J, Beale A.D. and Catchpole, M.A. (1996):
  A statistical algorithm for the early detection of outbreaks of
  infectious disease. J. R. Statist. Soc. A, 159, 547-563.

  Noufaily, A., Enki, D.G., Farrington, C.P., Garthwaite, P.,
  Andrews, N.J., Charlett, A. (2012): An improved algorithm for outbreak
  detection in multiple surveillance systems. Statistics in Medicine,
  32 (7), 1206-1222.

  Salmon, M., Schumacher, D., Stark, K., \enc{Höhle}{Hoehle}, M. (2015):
  Bayesian outbreak detection in the presence of reporting delays.
  Biometrical Journal, 57 (6), 1051-1067.
}
