\name{algo.glrnb}
\alias{algo.glrnb}
\encoding{latin1}

\title{Cound data regression charts}
\description{
 Count data regression charts for the monitoring of surveillance time series.
}
\usage{
algo.glrnb(disProgObj,control = list(range=range,c.ARL=5, 
         mu0=NULL, alpha=0, Mtilde=1, M=-1, change="intercept",
         theta=NULL,dir=c("inc","dec"),ret=c("cases","value")))
}
\arguments{
\item{disProgObj}{object of class \code{disProg} to do surveillance
    for}
\item{control}{A list controlling the behaviour of the algorithm
    \describe{
    \item{\code{range}}{vector of indices in the observed vector
        to monitor (should be consecutive)}
   
    \item{\code{mu0}}{A vector of in-control values of the mean of the negative binomial 
        distribution with the same length as \code{range}. If
        \code{NULL} the observed values in \code{1:(min(range)-1)} are
        used to estimate beta through a generalized linear model. To
        fine-tune the model one can instead specify \code{mu0} as a
        list with two components:
        \describe{
        \item{\code{S}}{number of harmonics to include}
        \item{\code{trend}}{include a term \code{t} in the GLM model}
        }}
    \item{\code{alpha}}{The (known) dispersion parameter of the
        negative binomial distribution. If \code{alpha=0} then the negative
        binomial distribution boils down to the Poisson distribution
        and a call of \code{algo.glrnb} is equivalent to a call to
        \code{algo.glrpois}. If \code{alpha=NULL} the parameter is calculated as part of the in-control estimation.}
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
        glr-statistic is computed (see below).}  
    } 
  } 
} 
\value{
\item{survRes}{ \code{algo.glrnb} returns a list of class
    \code{survRes} (surveillance result), which includes the alarm
    value for recognizing an outbreak (1 for alarm, 0 for no alarm),
    the threshold value for recognizing the alarm and the input object
    of class disProg. The \code{upperbound} slot of the object are
    filled with the current \eqn{GLR(n)} value or with the number of
    cases that are necessary to produce an alarm at any timpoint
    \eqn{<=n}. Both lead to the same alarm timepoints, but
    \code{"cases"} has an obvious interpretation.} 
}
% 
\details{ 
%
  This function implements the seasonal cound data chart based on
  generalized likelihood ratio (GLR) as described in the Hoehle and Paul
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
resetted starting from time \eqn{N+1}. Note that the same
\code{c.ARL} as before is used, but if \code{mu0} is different at
\eqn{N+1,N+2,\ldots} compared to time \eqn{1,2,\ldots} the run length
properties differ. Because \code{c.ARL} to obtain a specific ARL can
only be obtained my Monte Carlo simulation there is no good way to
update \code{c.ARL} automatically at the moment. Also, FIR GLR-detectors
might be worth considering.

At the moment, window limited ``\code{intercept}'' charts have not been
extensively tested and are at the moment not supported. As speed is
not an issue here this doesn't bother too much. Therefore, a value of
\code{M=-1} is always used in the intercept charts.

}
\seealso{
\code{\link{algo.rkiLatestTimepoint}}
}
\author{M. Hoehle}
\examples{
##Simulate data and apply the algorithm
S <- 1 ; t <- 1:120 ; m <- length(t)
beta <- c(1.5,0.6,0.6)
omega <- 2*pi/52
#log mu_{0,t}
alpha <- 0.2
base <- beta[1] + beta[2] * cos(omega*t) + beta[3] * sin(omega*t) 
#Generate example data with changepoint and tau=tau
tau <- 100
kappa <- 0.4
mu0 <- exp(base)
mu1 <- exp(base  + kappa) 

#Generate data
set.seed(42)
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
\source{
Count data regression charts for the monitoring of surveillance time
series (2008), M. \enc{Höhle}{Hoehle} and M. Paul, Computational Statistics and
Data Analysis, 52(9), pp. 4357--4368.

Poisson regression charts for the monitoring of surveillance time
series (2006), \enc{Höhle}{Hoehle}, M., SFB386 Discussion Paper 500.
}



% Hook

%Rnw file 
