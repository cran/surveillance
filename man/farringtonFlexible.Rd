\name{farringtonFlexible}
\alias{farringtonFlexible}
\encoding{latin1}

\title{Surveillance for Univariate Count Time Series Using an Improved Farrington Method}
\description{
%
  The function takes \code{range} values of the surveillance time
  series \code{sts} and for each time point uses a Poisson GLM with overdispersion to
  predict an upper bound on the number of counts according to the procedure by
  Farrington et al. (1996) and by Noufaily et al. (2012). This bound is then compared to the observed
  number of counts. If the observation is above the bound, then an alarm is raised.
  The implementation is illustrated in Salmon et al. (2016).
%
}
\usage{
farringtonFlexible(sts, control = list(
    range = NULL, b = 5, w = 3,
    reweight = TRUE, weightsThreshold = 2.58,
    verbose = FALSE, glmWarnings = TRUE,
    alpha = 0.05, trend = TRUE, pThresholdTrend = 0.05,
    limit54 = c(5,4), powertrans = "2/3",
    fitFun = "algo.farrington.fitGLM.flexible",
    populationOffset = FALSE,
    noPeriods = 1, pastWeeksNotIncluded = NULL,
    thresholdMethod = "delta"))
}
\arguments{
\item{sts}{object of class \code{\linkS4class{sts}} (including the \code{observed} and the \code{state} time series)}
\item{control}{Control object given as a \code{list} containing the following components:
    \describe{
    \item{\code{range}}{Specifies the index of all timepoints which
        should be tested. If range is \code{NULL} all possible timepoints are used.}
    \item{\code{b}}{How many years back in time to include when
        forming the base counts.}
    \item{\code{w}}{Window's half-size, i.e. number of weeks to include
        before and after the current week in each year.}
    \item{\code{reweight}}{Boolean specifying whether to perform reweighting step.}
     \item{\code{weightsThreshold}}{Defines the threshold for reweighting past outbreaks  using the Anscombe residuals
(1 in the original method, 2.58 advised in the improved method).}
      \item{\code{verbose}}{Boolean specifying whether to show extra debugging information.}
    \item{\code{glmWarnings}}{Boolean specifying whether to print warnings from the call to \code{glm}.}
    \item{\code{alpha}}{An approximate (one-sided) \eqn{(1-\alpha)\cdot 100\%}
        prediction interval is calculated unlike the original method where it was a two-sided interval. The upper limit of this interval
        i.e. the \eqn{(1-\alpha)\cdot 100\%} quantile serves as an upperbound.}
    \item{\code{trend}}{Boolean indicating whether a trend should be included and kept in
        case the conditions in the Farrington et. al. paper are met
        (see the results). If \code{false} then NO trend is fit.}
     \item{\code{pThresholdTrend}}{Threshold for deciding whether to keep trend in the model
(0.05 in the original method, 1 advised in the improved method).}
    \item{\code{limit54}}{Vector containing two numbers: \code{cases} and \code{period}. To avoid alarms in cases where the time series only
        has about almost no cases in the specific week the algorithm uses the following heuristic
        criterion (see Section 3.8 of the Farrington paper) to protect
        against low counts: no alarm is sounded if fewer than
        \eqn{\code{cases}=5} reports were received in the past \eqn{\code{period}=4}
        weeks. \code{limit54=c(cases,period)} is a vector allowing the
        user to change these numbers. Note: As of version 0.9-7 of the package the
        term "last" period of weeks includes the current week -
        otherwise no alarm is sounded for horrible large numbers if
        the four weeks before that are too low.}
    \item{\code{powertrans}}{Power transformation to apply to the
        data if the threshold is to be computed with  the method described
		in Farrington et al. (1996. Use either "2/3" for skewness correction (Default),
        "1/2" for variance stabilizing transformation or "none" for no
        transformation.}
    \item{\code{fitFun}}{String containing the name of the fit
        function to be used for fitting the GLM. The only current option is 
        "algo.farrington.fitGLM.flexible".}
    \item{\code{populationOffset}}{Boolean specifying whether to include 
                                 a population offset in the GLM. 
                                 The slot \code{sts@population} gives the population vector.}
    \item{\code{noPeriods}}{Number of levels in the factor allowing to use more baseline. If
    equal to 1 no factor variable is created, the set of reference values is defined as in
    Farrington et al (1996).}
    \item{\code{pastWeeksNotIncluded}}{Number of past weeks to ignore in the calculation.
      The default (\code{NULL}) means to use the value of \code{control$w}.
      Setting \code{pastWeeksNotIncluded=26} might be preferable
      (Noufaily et al., 2012).}
    \item{\code{thresholdMethod}}{Method to be used to derive the upperbound.
    Options are \code{"delta"} for the method described in Farrington et al. (1996),
    \code{"Noufaily"} for the method described in Noufaily et al. (2012),
    and \code{"muan"} for the method extended from Noufaily et al. (2012).}
      }
    }
}
\details{
  The following steps are performed according to the Farrington
  et al. (1996) paper.
\enumerate{
\item Fit of the initial model with intercept, time trend if \code{trend} is \code{TRUE},
seasonal factor variable if \code{noPeriod} is bigger than 1, and population offset if
\code{populationOffset} is \code{TRUE}. Initial estimation of mean and
  overdispersion.
\item Calculation of the weights omega (correction for past outbreaks) if \code{reweighting} is \code{TRUE}.
The threshold for reweighting is defined in \code{control}.
\item Refitting of the model
\item Revised estimation of overdispersion
\item Omission of the trend, if it is not significant
\item Repetition of the whole procedure
\item Calculation of the threshold value using the model to compute a quantile of the predictive distribution.
The method used depends on \code{thresholdMethod}, this can either be:  
\describe{
    \item{"delta"}{One assumes that the prediction error (or a transformation of the prediction
error, depending on \code{powertrans}), is normally distributed. The threshold is deduced from a quantile of
this normal distribution using the variance and estimate of the expected
count given by GLM, and the delta rule. The procedure takes into account both the estimation error (variance of the estimator
of the expected count in the GLM) and the prediction error (variance of the prediction error). This is the suggestion
in Farrington et al. (1996).}
    \item{"nbPlugin"}{One assumes that the new count follows
a negative binomial distribution parameterized by the expected count and the overdispersion
estimated in the GLM. The threshold is deduced from a quantile of this discrete distribution.
This process disregards the estimation error, though. This method was used in Noufaily, et al. (2012).}
    \item{"muan"}{One also uses the assumption of the negative
binomial sampling distribution but does not plug in the estimate of the expected count from the GLM,
instead one uses a quantile from the asymptotic normal distribution of the expected count estimated
in the GLM; in order to take into account both the estimation error and the prediction error.   }        }
\item Computation of exceedance score
}
Warning: monthly data containing the last day of each month as date should be analysed with \code{epochAsDate=FALSE} in the \code{sts} object. Otherwise February makes it impossible to find some reference time points. 
}
\value{
An object of class \code{sts} with the slots \code{upperbound} and \code{alarm} filled by appropriate output of the algorithm.
The \code{control} slot of the input \code{sts} is amended with the
following matrix elements, all with \code{length(range)} rows:
\describe{
\item{trend}{Booleans indicating whether a time trend was fitted for this time point.}
\item{trendVector}{coefficient of the time trend in the GLM for this time point. If no trend was fitted it is equal to NA.}
\item{pvalue}{probability of observing a value at least equal to the observation under the null hypothesis .}
\item{expected}{expectation of the predictive distribution for each timepoint. 
It is only reported if the conditions for raising an alarm are met (enough cases).}
\item{mu0Vector}{input for the negative binomial distribution to get the upperbound as a quantile 
(either a plug-in from the GLM or a quantile from the asymptotic normal distribution of the estimator)}
\item{phiVector}{overdispersion of the GLM at each timepoint.}
}
}

\keyword{classif}

\examples{

### DATA I/O ###
#Read Salmonella Agona data
data("salmonella.agona")

# Create the corresponding sts object from the old disProg object
salm <- disProg2sts(salmonella.agona)

### RUN THE ALGORITHMS WITH TWO DIFFERENT SETS OF OPTIONS ###           
# Farrington with old options
control1 <-  list(range=(260:312),
                  noPeriods=1,populationOffset=FALSE,
                  fitFun="algo.farrington.fitGLM.flexible",
                  b=4,w=3,weightsThreshold=1,
                  pastWeeksNotIncluded=3,
                  pThresholdTrend=0.05,trend=TRUE,
                  thresholdMethod="delta",alpha=0.1)
control2 <- list(range=(260:312),
                 noPeriods=10,populationOffset=FALSE,
                 fitFun="algo.farrington.fitGLM.flexible",
                 b=4,w=3,weightsThreshold=2.58,
                 pastWeeksNotIncluded=26,
                 pThresholdTrend=1,trend=TRUE,
                 thresholdMethod="delta",alpha=0.1)
salm1 <- farringtonFlexible(salm,control=control1)
salm2 <- farringtonFlexible(salm,control=control2)

### PLOT THE RESULTS ###
y.max <- max(upperbound(salm1),observed(salm1),upperbound(salm2),na.rm=TRUE)
 
plot(salm1, ylim=c(0,y.max), main='S. Newport in Germany', legend.opts=NULL)
lines(1:(nrow(salm1)+1)-0.5, 
      c(upperbound(salm1),upperbound(salm1)[nrow(salm1)]),
      type="s",col='tomato4',lwd=2)
lines(1:(nrow(salm2)+1)-0.5, 
      c(upperbound(salm2),upperbound(salm2)[nrow(salm2)]),
      type="s",col="blueviolet",lwd=2)

legend(0, 10, legend=c('Alarm','Upperbound with old options',
                       'Upperbound with new options'),
       pch=c(24,NA,NA),lty=c(NA,1,1),
       bg="white",lwd=c(2,2,2),col=c('red','tomato4',"blueviolet"))
}
\author{M. Salmon, M. \enc{Höhle}{Hoehle}}
\seealso{\code{\link{algo.farrington.fitGLM}},\code{\link{algo.farrington.threshold}}}
\keyword{classif}
\references{
  Farrington, C.P., Andrews, N.J, Beale A.D. and Catchpole, M.A. (1996):
  A statistical algorithm for the early detection of outbreaks of
  infectious disease. J. R. Statist. Soc. A, 159, 547-563.

  Noufaily, A., Enki, D.G., Farrington, C.P., Garthwaite, P., Andrews,
  N.J., Charlett, A. (2012): An improved algorithm for outbreak
  detection in multiple surveillance systems. Statistics in Medicine,
  32 (7), 1206-1222.

  Salmon, M., Schumacher, D. and \enc{Höhle}{Hoehle}, M. (2016):
  Monitoring count time series in \R: Aberration detection in public
  health surveillance. \emph{Journal of Statistical Software},
  \bold{70} (10), 1-35. \doi{10.18637/jss.v070.i10}
}
