\name{algo.glrpois}
\alias{algo.glrpois}
\title{Poisson regression charts}
\description{
 Poisson regression charts for the monitoring of surveillance time series
}
\usage{
algo.glrpois(disProgObj,control = list(range=range,c.ARL=5, S=1,
         mu0=NULL, Mtilde=1, M=-1, change="intercept",theta=NULL))
}
\arguments{
  \item{disProgObj}{object of class \code{disProg} to do surveillance
  for}
  \item{control}{A list controlling the behaviour of the algorithm.
                 \describe{
		 \item{range}{vector of indices in the observed vector
		 to monitor (should be consecutive)}
		 \item{S}{number of harmonics to include. Note:
		 Atm. only S=1 is supported.}
		 \item{beta}{A vector of regression coefficients to
		 use in the seasonal Poisson model. If \code{NULL} the
		 observed values in \code{1:min(range)} are used to
		 estimate beta through a generalized linear model.}
		 \item{c.ARL}{threshold in the GLR test,
		 i.e. \eqn{c_gamma}{c_{\gamma}}}
                 \item{Mtilde}{number of observations needed before we
		   have a full rank the typical setup for the
		   "intercept" and "epi" charts is Mtilde=1}
		 \item{M}{Number of time instances back in time in the
		 window-limited approach, i.e. the last value
		 considered is \eqn{\max{1,n-M}}. To always look back
		 until the first observation use \code{M=-1}.}
		 \item{change}{A string specifying the type of the
		 alternative. Currently the two choices are
		 \code{intercept} and \code{epi}. See the SFB
		 Discussion Paper 500 for details.}
		 \item{theta}{If \code{NULL} then the GLR scheme is
                 used. If not \code{NULL} the prespecified value for
                 \eqn{\kappa} or \eqn{lambda} is used in a recursive
                 LR scheme, which is faster.  }
		 }
  }
}
\value{
\item{survRes}{ \code{algo.prc} returns a list of class \code{survRes}
    (surveillance result), which includes the alarm value for
    recognizing an outbreak (1 for alarm, 0 for no alarm), the
    threshold value for recognizing the alarm and the input object of
    class disProg. The \code{upperbound} slot of the object are filled
    with the current \eqn{GLR(n)} value.  } 
}
%
\details{ 
%
  This function implements the seasonal Poisson charts based on
  generalized likelihood ratio (GLR) as described in the SFB
  Discussion Paper 500. A moving-window generalized likelihood ratio
  detector is used, i.e. the detector has the form
%
\deqn{N = inf(... >= c_gamma)}{N = \inf\left\{ n : \max_{1\leq k \leq
      n} \left[ \sum_{t=k}^n \log \left\{
        \frac{f_{\theta_1}(x_t|z_t)}{f_{\theta_0}(x_t|z_t)} \right\}
    \right] \geq c_\gamma \right\} }
%
where instead of \eqn{1<= k <= n}{1\leq k \leq n} the GLR statistic is computed
for all \eqn{k \in \{n-M, \ldots, n-\tilde{M}+1\}}. To achieve the
typical behaviour from \eqn{1\leq k\leq n} use \code{Mtilde=1} and \code{M=-1}.

At the moment, window limited ``intercept'' charts have not been
tested and are atm not supported. As speed is not an issue here this
doesn't bother too much. Therefore, a value of \code{M=-1} is always
in the intercept charts.

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
base <- beta[1] + beta[2] * cos(omega*t) + beta[3] * sin(omega*t) 
#Generate example data with changepoint and tau=tau
tau <- 100
kappa <- 0.4
mu0 <- exp(base)
mu1 <- exp(base  + kappa) 

#Generate data
set.seed(42)
x <- rpois(length(t),mu0*(exp(kappa)^(t>=tau)))
s.ts <- create.disProg(week=1:length(t),observed=x,state=(t>=tau))

#Plot the data
plot(s.ts,legend=NULL,xaxis.years=FALSE)

#Run 
cntrl = list(range=t,c.ARL=5, S=1,Mtilde=1, mu0=mu0,change="intercept")
glr.ts <- algo.glrpois(s.ts,control=cntrl)
lr.ts  <- algo.glrpois(s.ts,control=c(cntrl,theta=0.4))

plot(glr.ts,xaxis.years=FALSE)

}

\keyword{}
\source{
Poisson regression charts for the monitoring of surveillance time
series (2006), Höhle, M., SFB386 Discussion Paper 500.
}

