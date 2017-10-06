\name{algo.glrnb}
\alias{algo.glrnb}
\alias{algo.glrpois}
\encoding{latin1}

\title{Count Data Regression Charts}

\description{
 Count data regression charts for the monitoring of surveillance time
 series as proposed by \enc{Höhle}{Hoehle} and Paul (2008).
 The implementation is described in Salmon et al. (2016).
}

\usage{
algo.glrnb(disProgObj, control = list(range=range, c.ARL=5,
           mu0=NULL, alpha=0, Mtilde=1, M=-1, change="intercept",
           theta=NULL, dir=c("inc","dec"),
           ret=c("cases","value"), xMax=1e4))

algo.glrpois(disProgObj, control = list(range=range, c.ARL=5,
             mu0=NULL, Mtilde=1, M=-1, change="intercept",
             theta=NULL, dir=c("inc","dec"),
             ret=c("cases","value"), xMax=1e4))
}

\arguments{
\item{disProgObj}{object of class \code{disProg} to do surveillance
    for}
\item{control}{A list controlling the behaviour of the algorithm
    \describe{
    \item{\code{range}}{vector of indices in the observed vector
        to monitor (should be consecutive)}

    \item{\code{mu0}}{A vector of in-control values of the mean of the
        Poisson / negative binomial
        distribution with the same length as \code{range}. If
        \code{NULL} the observed values in \code{1:(min(range)-1)} are
        used to estimate the beta vector through a generalized linear
        model. To
        fine-tune the model one can instead specify \code{mu0} as a
        list with two components:
        \describe{
          \item{\code{S}}{integer number of harmonics to include
        (typically 1 or 2)}
          \item{\code{trend}}{A Boolean indicating whether to include a term \code{t} in the GLM model}
	}
	The fitting is controlled by the \code{estimateGLRNbHook}
        function. The in-control mean model is re-fitted after every
	alarm. The fitted models can be found as a list \code{mod} in
        the \code{control} slot after the call.

	Note: If a value for \code{alpha} is given, then the inverse of
        this value is used as fixed \code{theta} in a
        \code{\link[MASS]{negative.binomial}} \code{glm}.
        If \code{is.null(alpha)} then the parameter is
        estimated as well (using \code{\link[MASS]{glm.nb}}) --
        see the description of this parameter for details.
      }
    \item{\code{alpha}}{The (known) dispersion parameter of the negative
        binomial distribution, i.e. the parametrization of the negative
        binomial is such that the variance is \eqn{mean +
        alpha*mean^2}{mean + \alpha*mean^2}. Note: This parametrization
        is the inverse of the shape parametrization used in R -- for
        example in \code{dnbinom} and \code{glr.nb}. Hence, if
        \code{alpha=0} then the negative binomial distribution boils
        down to the Poisson distribution and a call of \code{algo.glrnb}
        is equivalent to a call to \code{algo.glrpois}. If
        \code{alpha=NULL} the parameter is calculated as part of the
        in-control estimation. However, the parameter is estimated only
        once from the first fit. Subsequent fittings are only for the
        parameters of the linear predictor with \code{alpha} fixed.}
    \item{\code{c.ARL}}{threshold in the GLR test,
        i.e. \eqn{c_{\gamma}}{c_gamma}}
    \item{\code{Mtilde}}{number of observations needed before we
        have a full rank the typical setup for the
        "\code{intercept}" and "\code{epi}" charts is \code{Mtilde=1}}
    \item{\code{M}}{number of time instances back in time in the
        window-limited approach, i.e. the last value
        considered is \eqn{\max{1,n-M}}. To always look back
        until the first observation use \code{M=-1}.}
    \item{\code{change}}{a string specifying the type of the
        alternative. Currently the two choices are
        \code{intercept} and \code{epi}. See the SFB
        Discussion Paper 500 for details.}
    \item{\code{theta}}{if \code{NULL} then the GLR scheme is
        used. If not \code{NULL} the prespecified value for
        \eqn{\kappa} or \eqn{\lambda} is used in a recursive
        LR scheme, which is faster.  }
    \item{\code{dir}}{a string specifying the direction of testing in
        GLR scheme. With \code{"inc"} only increases in \eqn{x} are
        considered in the GLR-statistic, with \code{"dec"} decreases
        are regarded. }
    \item{\code{ret}}{a string specifying the type of
        \code{upperbound}-statistic that is returned. With
        \code{"cases"} the number of cases that would have been
        necessary to produce an alarm or with \code{"value"} the
        GLR-statistic is computed (see below).}
      \item{\code{xMax}}{Maximum value to try for x to see if this is
        the upperbound number of cases before sounding an alarm
        (Default: 1e4). This only applies for the GLR using the NegBin
        when \code{ret="cases"} -- see details.}
    }
  }
}

\value{
  \code{algo.glrpois} simply calls \code{algo.glrnb} with
  \code{control$alpha} set to 0.

  \code{algo.glrnb} returns a list of class
    \code{survRes} (surveillance result), which includes the alarm
    value for recognizing an outbreak (1 for alarm, 0 for no alarm),
    the threshold value for recognizing the alarm and the input object
    of class disProg. The \code{upperbound} slot of the object are
    filled with the current \eqn{GLR(n)} value or with the number of
    cases that are necessary to produce an alarm at any time point
    \eqn{<=n}. Both lead to the same alarm timepoints, but
    \code{"cases"} has an obvious interpretation.
}

\details{

  This function implements the seasonal count data chart based on
  generalized likelihood ratio (GLR) as described in the \enc{Höhle}{Hoehle} and Paul
  (2008) paper. A moving-window generalized likelihood ratio
  detector is used, i.e. the detector has the form
%
\deqn{N = \inf\left\{ n : \max_{1\leq k \leq
      n} \left[ \sum_{t=k}^n \log \left\{
        \frac{f_{\theta_1}(x_t|z_t)}{f_{\theta_0}(x_t|z_t)} \right\}
    \right] \geq c_\gamma \right\} }{N = inf(... >= c_gamma)}
%
where instead of \eqn{1\leq k \leq n}{1<= k <= n} the GLR statistic is
computed for all \eqn{k \in \{n-M, \ldots, n-\tilde{M}+1\}}{k \in \{n-M, \ldots, n-Mtilde+1\}}. To
achieve the typical behaviour from \eqn{1\leq k\leq n}{1<= k <= n} use
\code{Mtilde=1} and \code{M=-1}.

So \eqn{N} is the time point where the GLR statistic is above the
threshold the first time: An alarm is given and the surveillance is
reset starting from time \eqn{N+1}. Note that the same
\code{c.ARL} as before is used, but if \code{mu0} is different at
\eqn{N+1,N+2,\ldots} compared to time \eqn{1,2,\ldots} the run length
properties differ. Because \code{c.ARL} to obtain a specific ARL can
only be obtained my Monte Carlo simulation there is no good way to
update \code{c.ARL} automatically at the moment. Also, FIR GLR-detectors
might be worth considering.

In case \code{is.null(theta)} and \code{alpha>0} as well as
\code{ret="cases"} then a brute-force search is conducted for each time
point in range in order to determine the number of cases necessary
before an alarm is sounded. In case no alarm was sounded so far by time
\eqn{t}, the function increases \eqn{x[t]} until an alarm is sounded any
time before time point \eqn{t}. If no alarm is sounded by \code{xMax}, a return value
of 1e99 is given. Similarly, if an alarm was sounded by time \eqn{t} the
function counts down instead. Note: This is slow experimental code!

At the moment, window limited ``\code{intercept}'' charts have not been
extensively tested and are at the moment not supported. As speed is
not an issue here this doesn't bother too much. Therefore, a value of
\code{M=-1} is always used in the intercept charts.

}

\author{M. \enc{Höhle}{Hoehle} with contributions by V. Wimmer}

\examples{
##Simulate data and apply the algorithm
S <- 1 ; t <- 1:120 ; m <- length(t)
beta <- c(1.5,0.6,0.6)
omega <- 2*pi/52
#log mu_{0,t}
base <- beta[1] + beta[2] * cos(omega*t) + beta[3] * sin(omega*t)
#Generate example data with changepoint and tau=tau
tau <- 100
kappa <- 0.4
mu0 <- exp(base)
mu1 <- exp(base  + kappa)


## Poisson example
#Generate data
set.seed(42)
x <- rpois(length(t),mu0*(exp(kappa)^(t>=tau)))
s.ts <- create.disProg(week=1:length(t),observed=x,state=(t>=tau))
#Plot the data
plot(s.ts,legend=NULL,xaxis.years=FALSE)
#Run
cntrl = list(range=t,c.ARL=5, Mtilde=1, mu0=mu0,
             change="intercept",ret="value",dir="inc")
glr.ts <- algo.glrpois(s.ts,control=cntrl)
plot(glr.ts,xaxis.years=FALSE)
lr.ts  <- algo.glrpois(s.ts,control=c(cntrl,theta=0.4))
plot(lr.ts,xaxis.years=FALSE)


## NegBin example
#Generate data
set.seed(42)
alpha <- 0.2
x <- rnbinom(length(t),mu=mu0*(exp(kappa)^(t>=tau)),size=1/alpha)
s.ts <- create.disProg(week=1:length(t),observed=x,state=(t>=tau))

#Plot the data
plot(s.ts,legend=NULL,xaxis.years=FALSE)

#Run GLR based detection
cntrl = list(range=t,c.ARL=5, Mtilde=1, mu0=mu0, alpha=alpha,
             change="intercept",ret="value",dir="inc")
glr.ts <- algo.glrnb(s.ts,control=c(cntrl))
plot(glr.ts,xaxis.years=FALSE)

#CUSUM LR detection with backcalculated number of cases
cntrl2 = list(range=t,c.ARL=5, Mtilde=1, mu0=mu0, alpha=alpha,
              change="intercept",ret="cases",dir="inc",theta=1.2)
glr.ts2 <- algo.glrnb(s.ts,control=c(cntrl2))
plot(glr.ts2,xaxis.years=FALSE)

}

\keyword{classif}

\references{
  \enc{Höhle}{Hoehle}, M. and Paul, M. (2008):
  Count data regression charts for the monitoring of surveillance time
  series. Computational Statistics and Data Analysis, 52 (9), 4357-4368.

  Salmon, M., Schumacher, D. and \enc{Höhle}{Hoehle}, M. (2016):
  Monitoring count time series in \R: Aberration detection in public
  health surveillance. \emph{Journal of Statistical Software},
  \bold{70} (10), 1-35. \doi{10.18637/jss.v070.i10}
}
