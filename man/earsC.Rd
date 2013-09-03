\name{earsC}
\alias{earsC}
\encoding{latin1}

\title{Surveillance for a count data time series using the EARS C1, C2 or C3 method.}
\description{
%
  The function takes \code{range} values of the surveillance time
  series \code{sts} and for each time point computes a threshold for the number of counts
   based on values from the recent past.
  This is then compared to the observed
  number of counts. If the observation is above a specific quantile of
  the prediction interval, then an alarm is raised.  This method is especially useful
  for data without many reference values, since it only needs counts from the recent past.
%
}
\usage{ earsC(sts, control = list(range = NULL, method = "C1",
                            alpha = 0.001))
}
\arguments{
\item{sts}{object of class sts (including the \code{observed} and the \code{state} time series) , which is to be monitored.}
\item{control}{Control object
    \describe{
    \item{\code{range}}{Specifies the index of all timepoints which
        should be tested. If \code{range} is \code{NULL} the maximum number
        of possible timepoints is used. This number depends on the method chosen.
        For C1 all timepoints from timepoint 8 can be assessed, for C2 from 
        timepoint 10 and for C3 from timepoint 12.}
    \item{\code{method}}{String indicating which method to use: \cr 
    "C1" for EARS C1-MILD method, "C2" for EARS C2-MEDIUM method,
    "C3" for EARS C3-HIGH method. By default if \code{method} is  \code{NULL} C1 is chosen.}

    \item{\code{alpha}}{An approximate (two-sided) \eqn{(1-\alpha)\cdot 100\%} prediction
     interval is calculated. By default if \code{alpha} is \code{NULL} 0.001 is assumed
     for C1 and C2 whereas 0.025 is assumed for C3. These different choices are the one made at the CDC.}
        %
      }
    }
}
\details{
  The three methods are different in terms of baseline used for calculation of 
  the expected value and in terms of method for calculating the expected value:
 \itemize{
 \item in C1 and C2 the expected value is the moving average of counts over 
 the sliding window of the baseline and the prediction interval depends on the 
 standard derivation of counts over this window. They can be considered as
 Shewhart control charts with a small sample used for calculations.
 \item in C3 the expected value is based on the sum over 3 timepoints 
 (assessed timepoints and the two previous timepoints) of the discrepancy 
 between observations and predictions, predictions being calculated with C2 method.
This method shares a common point with CUSUM method 
(adding discrepancies between predictions and observations over several timepoints) 
but is not a CUSUM  (sum over 3 timepoints, not accumulation over a whole range),
 even if it sometimes presented as such.
}
Here is what the function does for each method:
\enumerate{
\item For C1 the baseline are the 7 timepoints before the assessed timepoint t,
 t-7 to t-1. The expected value is the mean of the baseline. An approximate
  (two-sided) \eqn{(1-\alpha)\cdot 100\%} prediction interval is calculated based on the
   assumption that  the difference between the expected value and the observed 
   value divided by the standard derivation of counts over the sliding window,
    called \eqn{C_1(t)}, follows a standard normal distribution in the absence
    of outbreaks:
\deqn{C_1(t)= \frac{Y(t)-\bar{Y}_1(t)}{S_1(t)},}
where
\deqn{\bar{Y}_1(t)= \frac{1}{7} \sum_{i=t-1}^{t-7} Y(i)}
and
\deqn{ S^2_1(t)= \frac{1}{6} \sum_{i=t-1}^{t-7} [Y(i) - \bar{Y}_1(i)]^2.}
Then under the null hypothesis of no outbreak,
 \deqn{C_1(t) \mathcal \sim {N}(0,1)}
An alarm is raised if \deqn{C_1(t)\ge z_{1-\alpha}}
with \eqn{z_{1-\alpha}} the \eqn{(1-\alpha)^{th}} quantile of the centered
reduced normal law. \cr

The upperbound \eqn{U_1(t)} is then defined by:
 \deqn{U_1(t)= \bar{Y}_1(t) + z_{1-\alpha}S_1(t).}
\item C2 is very close to C1 apart from a 2-day lag in the baseline definition.
Indeed for C2 the baseline are 7 timepoints with a 2-day lag before the assessed
timepoint t, t-9 to t-3. The expected value is the mean of the baseline. An approximate
  (two-sided) \eqn{(1-\alpha)\cdot 100\%} prediction interval is calculated based on the
   assumption that  the difference between the expected value and the observed 
   value divided by the standard derivation of counts over the sliding window,
    called \eqn{C_2(t)}, follows a standard normal distribution in the absence
    of outbreaks:
\deqn{C_2(t)= \frac{Y(t)-\bar{Y}_2(t)}{S_2(t)},} 
where
\deqn{\bar{Y}_2(t)= \frac{1}{7} \sum_{i=t-3}^{t-9} Y(i)}
and
\deqn{ S^2_2(t)= \frac{1}{6} \sum_{i=t-3}^{t-9} [Y(i) - \bar{Y}_2(i)]^2.}
Then under the null hypothesis of no outbreak,
 \deqn{C_2(t) \mathcal \sim {N}(0,1)}
An alarm is raised if \deqn{C_2(t)\ge z_{1-\alpha},}
with \eqn{z_{1-\alpha}} the \eqn{(1-\alpha)^{th}} quantile of the centered
reduced normal law. \cr

The upperbound \eqn{U_2(t)} is then defined by:
 \deqn{U_2(t)= \bar{Y}_2(t) + z_{1-\alpha}S_2(t).}
 
 \item C3 is quite different from the two other methods but it is based on C2.
 Indeed it uses \eqn{C_2(t)} from timepoint t and the two previous timepoints.
 This means the baseline are timepoints t-11 to t-3.
 The statistic \eqn{C_3(t)} is the sum of discrepancies between observations and 
 predictions.
\deqn{C_3(t)= \sum_{i=t}^{t-2} \max(0,C_2(i)-1)} 
 Then under the null hypothesis of no outbreak,
 \deqn{C_3(t) \mathcal \sim {N}(0,1)}
An alarm is raised if \deqn{C_3(t)\ge z_{1-\alpha},}
with \eqn{z_{1-\alpha}} the \eqn{(1-\alpha)^{th}} quantile of the centered
reduced normal law. \cr

The upperbound \eqn{U_3(t)} is then defined by:
 \deqn{U_3(t)= \bar{Y}_2(t) + S_2(t)\left(z_{1-\alpha}-\sum_{i=t-1}^{t-2} \max(0,C_2(i)-1)\right).}

}
}
\value{
 An object of class \code{sts} with the slots \code{upperbound} and \code{alarm} filled
 by the chosen method.
}
\examples{

#Sim data and convert to sts object
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208, A = 1,
                              alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.7)
stsObj = disProg2sts( disProgObj)


#Call function and show result
res1 <- earsC(stsObj, control = list(range = 20:208,method="C1"))
plot(res1,legend.opts=list(horiz=TRUE,x="topright"),dx.upperbound=0)


# compare upperbounds depending on alpha
res3 <- earsC(stsObj, control = list(range = 20:208,method="C3",alpha = 0.001))
plot(res3@upperbound,t='l')
res3 <- earsC(stsObj, control = list(range = 20:208,method="C3"))
lines(res3@upperbound,col='red')



}

\author{M. Salmon}
\keyword{classif}
\source{
  Fricker, R.D., Hegler, B.L, and Dunfee, D.A. (2008). Comparing syndromic surveillance detection methods: EARS versus a CUSUM-based methodology,
  27:3407-3429, Statistics in medicine.
}

