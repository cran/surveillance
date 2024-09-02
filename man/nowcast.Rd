\encoding{latin1}
\name{nowcast}
\alias{nowcast}
%Internal functions
%\alias{dist.median}
%\alias{outside.ci}
%\alias{logS}
%\alias{RPS}

\title{
  Adjust a univariate time series of counts for observed
  but-not-yet-reported events
}

\description{
  Nowcasting can help to obtain up-to-date information on trends during
  a situation where reports about events arrive with delay. For example
  in public health reporting, reports about important
  indicators (such as occurrence of cases) are prone to be delayed due to
  for example manual quality checking and reporting system
  hierarchies. Altogether, the delays are subject to a delay distribution,
  which may, or may not, vary over time.
}

\usage{
nowcast(now, when, data, dEventCol="dHospital", dReportCol="dReport",
        method=c("bayes.notrunc", "bayes.notrunc.bnb", "lawless",
                 "bayes.trunc", "unif", "bayes.trunc.ddcp"),
        aggregate.by="1 day",
        D=15,
        m=NULL, m.interpretation=c("hoehle_anderheiden2014", "lawless1994"),
        control=list(
            dRange=NULL, alpha=0.05, nSamples=1e3,
            N.tInf.prior=c("poisgamma","pois","unif"),
            N.tInf.max=300, gd.prior.kappa=0.1,
            ddcp=list(ddChangepoint=NULL,
                      cp_order=c("zero","one"),
                      Wextra=NULL,
                      logLambda=c("iidLogGa","tps","rw1","rw2"),
                      responseDistr=c("poisson", "negbin"),
                      mcmc=c(burnin=2500, sample=10000, thin=1, adapt=1000,
                             store.samples=FALSE)),
            score=FALSE, predPMF=FALSE))
}

\arguments{
  \item{now}{
    an object of class \code{Date} denoting the day at which to do the
    nowcast. This corresponds to \eqn{T} in the notation of
    \enc{Höhle}{Hoehle} and an der Heiden (2014).
  }
  \item{when}{a vector of \code{Date} objects denoting the day(s) for which
    the projections are to be done. One needs to ensure that each
    element in \code{when} is smaller or equal to \code{now}.
  }
  \item{data}{A data frame with one row per case -- for each case on needs
    information on the day of the event (e.g. hospitalization) and the
    day of report of this event.
  }
  \item{dEventCol}{The name of the column in \code{data} which contains the
    date of the event, e.g. hospitalization. Default: \code{"dHospital"}.
  }
  \item{dReportCol}{Name of the column in \code{data} containing the date at
    which the report arrives at the respective register. Default:
    \code{"dReport"}.
  }
  \item{method}{A vector of strings denoting the different methods for doing
    the nowcasting. Note that results of the first name in this list are
    officially
    returned by the function. However, it is possible to specify several
    methods here, e.g., in order to compare score evaluations. Details of
    the methods are described in \enc{Höhle}{Hoehle} and an der Heiden (2014).
    \describe{
      \item{\code{"unif"}}{}
      \item{\code{"bayes.notrunc"}}{A Bayesian procedure ignoring
	truncation.}
      \item{\code{"bayes.notrunc.bnb"}}{A fast Bayesian procedure ignoring
	truncation and which calculates the adjustment per-time
	(i.e. ignoring other delays) using the negative binomial.}
      \item{\code{"lawless"}}{A discretized version of the Gaussian
	predictive distribution suggested in Lawless (1994).}
      \item{\code{"bayes.trunc"}}{Bayesian method based on the generalized
	Dirichlet distribution, which is the conjugate prior-posterior for the
	delay distribution PMF under right-truncated sampling as shown in HadH
	(2014).}
      \item{\code{"bayes.trunc.ddcp"}}{Fully Bayesian method allowing for
	change-points in the delay distribution, e.g., due to speed-ups
	in the reporting process. A discrete-survival
	model is used for the delay distribution. Details of
	the methods are described in HadH (2014). Note: This method
	requires that the JAGS program is installed on the system.}
    }
  }
  \item{aggregate.by}{Time scale used for the temporal aggregation of
    the records in the data \code{data}. See
    \code{\link{linelist2sts}} and \code{\link{seq.Date}} for further
  information.}
  \item{D}{Maximum possible or maximum relevant delay (unit:
    \code{aggregate.by}). Default: 15.}
  \item{m}{Size of the moving window for the estimation of the delay
    distribution. Default: \code{NULL}, i.e. take all values at all
  times. Otherwise: a positive integer equal to or greater than \code{D}
  such that only values from a sliding window are used. The shape of the
  window depends on the value of \code{m.interpretation}.}
  \item{m.interpretation}{This parameter controls the interpretation of
  the sliding window used to estimate the delay distribution. If
  \code{m.interpretation="hoehle_anderheiden2014"} (Default) then the
  sliding window is defined as a horizontal cut in the reporting
  triangle, i.e. the values for the delay estimation originate from
  reports occurring during \code{(now-m):now}. This means that the estimation of long delays is
  based on fewer observations than the estimation of the short delays,
  hence, the long delay estimates are subject to more variability. If
  for example \eqn{m=D} then the estimate for a delay of \eqn{d=D} is
  based on only one observation.
  The advantage of this choice is that one explicitly knows which time
  period all observations originate from. For details see Section 3 of
  \enc{Höhle}{Hoehle} and an der Heiden (2014).    

  Alternatively, when \code{m.interpretation}="lawless1994", the cut in
  the reporting triangle is made such that each delay \code{d} is
  estimated based on the same number of observations (\eqn{m+1}). This
  means that in order to estimate the delay for \eqn{d} days, a
  sliding rectangle of length \eqn{m+1} containing the reports which
  occurred during \code{(now-m-d):now}. See Fig. 2 in Lawless (1994) for
  details. Note: A warning is given is
  \code{method="lawless"}, but \code{m.interpretation} is not.}
  \item{control}{A list with named arguments controlling the
    functionality of the nowcasting.
    \describe{
      \item{dRange}{Default: \code{NULL}. In this case the
	\code{dEventCol} column is used to extract the first and last
	available in \code{data}.}
      \item{alpha}{Equal tailed (1-\eqn{\alpha}{alpha})*100\% prediction
	intervals are calculated. Default: 0.05.}
      \item{nSamples}{Number of PMF samples in the \code{bayes.*}
	procedures. Note: Entire vectors containing the PMF on the grid from 0
	to \code{N.tInf.max} are drawn and which are then combined. The
	argument does not apply to the \code{bayes.trunc.ddcp} method.}
      \item{N.tInf.prior}{Prior distribution of
	\eqn{N(t,\infty)}{N(t,Inf)}. Applies only to the \code{bayes.*} except
	\code{bayes.bayes.ddcp} methods. See example on how to control the
	distribution parameters.}
      \item{N.tInf.max}{Limit of the support of \eqn{N(t,\infty)}{N(t,Inf)}. The
	value needs to be high enough such that at this limit only little of
	the predictive distribution is right-truncated. Default: 300.}
      \item{gd.prior.kappa}{Concentration parameter for the Dirichlet
	prior for the delay distribution on \eqn{0,...,D}. Default: 0.1. Note:
	The procedure is quite sensitive to this parameter in case only few
	cases are available.}
      \item{ddcp}{A list specifying the change point model for the delay
	distribution. This method should only be used if detailed information
	about changes in the delay distribution are available as, e.g., in the
	case of the STEC O104:H4 outbreak. The components are as
	follows:
	\describe{
	  \item{\code{ddChangepoint}}{Vector of Date objects corresponding
	    to the changepoints}
	  \item{\code{cp_order}}{Either \code{"zero"} (Default) or
  \code{"one"}. This is the degree of the TPS spline for the baseline
  hazard, which is formed by the changepoints. Order zero corresponds to
  the dummy variables of the change-points being simply zero or one. In
  case a 1st order polynomial is chosen, this allows the delay
  distribution to change towards faster or slow reporting as time
  progresses (until the next change-point). The later can be helpful in
  very dynamic epidemic situations where a lot of cases suddenly appear
  overwhelming the surveillance system infrastructure.}
	  \item{\code{Wextra}}{An additional design matrix part to be joined
  onto the part originating from the change-points. Altogether, the
  column bind of these two quantities will be \eqn{W_{t,d}}. This allows
  one to include, e.g., day of the week effects or holidays.}
	  \item{\code{logLambda}}{Prior on the spline. One of
	    \code{c("iidLogGa","tps","rw1","rw2")}.}
	  \item{\code{respDistr}}{Response distribution of \eqn{n_{t,d}}
  in the reporting triangle. Default is \code{"poisson"}. An
  experimental alternative is to use \code{"negbin"}.}
	  \item{\code{tau.gamma}}{}
	  \item{\code{eta.mu}}{Vector of coefficients describing the
  mean of the prior normal distribution of the regression effects in the
  discrete time survival model.}
	  \item{\code{eta.prec}}{A precision matrix for the regression
  effects in the discrete time survival model.}
	  \item{\code{mcmc}}{A named vector of length 5 containing
	    burn-in (default: 2500), number of samples (10000),
            thinning (1) and adaptation (1000) for the
            three MCMC chains which are ran. The values are passed on to
  \code{\link[runjags]{run.jags}}. The fifth argument \code{store.samples}
  denotes if the output of the JAGS sampling should be included as part of the
  returned \code{stsNC} object. Warning: If \code{TRUE} (Default: \code{FALSE})
  the size of the returned object might increase substantially.}
	}
      }
      \item{score}{Compute scoring rules. Default: \code{FALSE}. The
	computed scores are found in the \code{SR} slot of the result.}
      \item{predPMF}{Boolean whether to return the probability mass
	functions of the individual forecasts (Default: \code{FALSE}). The
	result can be found in the \code{control} slot of the return object.}
    }
  }
}

\details{
  The methodological details of the nowcasting procedures are described in
  \enc{Höhle}{Hoehle} M and an der Heiden M (2014).
}

\value{
  \code{nowcast} returns an object of \code{"\linkS4class{stsNC}"}. The
  \code{upperbound} slot contains the median of the method specified at
  the first position the argument \code{method}. The slot \code{pi} (for
  prediction interval)
  contains the equal tailed (1-\eqn{\alpha}{alpha})*100\% prediction
  intervals, which are calculated based on the predictive distributions
  in slot \code{predPMF}.
  Furthermore, slot \code{truth} contains an \code{sts} object
  containing the true number of cases (if possible to compute it is based on
  the data in \code{data}). Finally, slot \code{SR} contains the results
  for the proper scoring rules (requires truth to be calculable).
}

\references{
  \enc{Höhle}{Hoehle}, M. and an der Heiden, M. (2014): Bayesian nowcasting
  during the STEC O104:H4 outbreak in Germany, 2011. \emph{Biometrics}
  70(4):993-1002. \doi{10.1111/biom.12194}.\cr
  A preprint is available as
  \url{https://staff.math.su.se/hoehle/pubs/hoehle_anderheiden2014-preprint.pdf}.

  \enc{Günther}{Guenther}, F. and Bender, A. and Katz, K. and
  \enc{Küchenhoff}{Kuechenhoff}, H. and \enc{Höhle}{Hoehle}, M. (2020):
  Nowcasting the COVID-19 pandemic in Bavaria.
  \emph{Biometrical Journal}. \doi{10.1002/bimj.202000112}\cr
  Preprint available at \doi{10.1101/2020.06.26.20140210}.
}


\author{
  Michael \enc{Höhle}{Hoehle}
}

\note{
  Note: The \code{bayes.trunc.ddcp} uses the JAGS software together with
  the \R package \pkg{runjags} to handle the parallelization of
  the MCMC using the \code{"rjparallel"} method of
  \code{\link[runjags]{run.jags}}, which additionally requires the
  \pkg{rjags} package. You need to manually install
  JAGS on your computer for the package to work -- see
  \url{https://mcmc-jags.sourceforge.io/}
  and the documentation of \pkg{runjags} for details.

  Note: The function is still under development and might change in the
  future. Unfortunately, little emphasis has so far been put on making
  the function easy to understand and use.
}

\examples{
data("husO104Hosp")

#Extract the reporting triangle at a specific day
t.repTriangle <- as.Date("2011-07-04")

#Use 'void' nowcasting procedure (we just want the reporting triangle)
nc <- nowcast(now=t.repTriangle,when=t.repTriangle,
              dEventCol="dHosp",dReportCol="dReport",data=husO104Hosp,
              D=15,method="unif")

#Show reporting triangle
reportingTriangle(nc)

#Perform Bayesian nowcasting assuming the delay distribution is stable over time
nc.control <- list(N.tInf.prior=structure("poisgamma",
                                mean.lambda=50,var.lambda=3000),
                                nSamples=1e2)

t.repTriangle <- as.Date("2011-06-10")
when <- seq(t.repTriangle-3,length.out=10,by="-1 day")
nc <- nowcast(now=t.repTriangle,when=when,
              dEventCol="dHosp",dReportCol="dReport",data=husO104Hosp,
              D=15,method="bayes.trunc",control=nc.control)

#Show time series and posterior median forecast/nowcast
plot(nc,xaxis.tickFreq=list("\%d"=atChange,"\%m"=atChange),
     xaxis.labelFreq=list("\%d"=at2ndChange),xaxis.labelFormat="\%d-\%b",
     xlab="Time (days)",lty=c(1,1,1,1),lwd=c(1,1,2))

\dontrun{
### Using runjags to do a Bayesian model with changepoint(s)
### -- this might take a while
nc.control.ddcp <- modifyList(nc.control,
                    list(gd.prior.kappa=0.1,
                         ddcp=list(ddChangepoint=as.Date(c("2011-05-23")),
                             logLambda="tps",
                             tau.gamma=1,
                             mcmc=c(burnin=1000,sample=1000,thin=1,
                                    adapt=1000,store.samples=FALSE))))

nc.ddcp <- nowcast(now=t.repTriangle,when=when,
               dEventCol="dHosp",dReportCol="dReport",
               data=husO104Hosp, aggregate.by="1 day",
               method="bayes.trunc.ddcp", D=15,
                   control=nc.control.ddcp)

plot(nc.ddcp,legend.opts=NULL,
     xaxis.tickFreq=list("\%d"=atChange,"\%m"=atChange),
     xaxis.labelFreq=list("\%d"=at2ndChange),xaxis.labelFormat="\%d-\%b",
     xlab="Time (days)",lty=c(1,1,1,1),lwd=c(1,1,2))

lambda <- attr(delayCDF(nc.ddcp)[["bayes.trunc.ddcp"]],"model")$lambda
showIdx <- seq(which( max(when) == epoch(nc.ddcp))) #seq(ncol(lambda))
matlines( showIdx,t(lambda)[showIdx,],col="gray",lwd=c(1,2,1),lty=c(2,1,2))
legend(x="topright",c(expression(lambda(t)),"95\% CI"),col="gray",lwd=c(2,1),lty=c(1,2))
}
}

\keyword{models}
