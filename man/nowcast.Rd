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
  which may or may not vary over time.
}

\usage{
nowcast(now,when,data,dEventCol="dHospital",dReportCol="dReport",
        method=c("bayes.notrunc","bayes.notrunc.bnb","lawless","bayes.trunc",
                 "unif","bayes.trunc.ddcp"),
        aggregate.by="1 day",
        D=15, m=NULL,
        control=list(
            dRange=NULL,alpha=0.05,nSamples=1e3,
            N.tInf.prior=c("poisgamma","pois","unif"),
            N.tInf.max=300, gd.prior.kappa=0.1,
            ddcp=list(ddChangepoint=NULL,
                logLambda=c("iidLogGa","tps","rw1","rw2"),
                tau.gamma=1,eta.mu=NULL, eta.prec=NULL,
                mcmc=c(burnin=2500,sample=10000,thin=1)),
            score=FALSE,predPMF=FALSE))
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
  times. Otherwise: an integer such that values from \code{(now-m):now}
  are used.}
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
	The procedure is quite sensitive to this parameters in case only few
	cases are available.}
      \item{ddcp}{A list specifying the change point model for the delay
	distribution. This method should only be used if detailed information
	about changes in the delay distribution are available as, e.g., in the
	case of the STEC O104:H4 outbreak. The components are as
	follows:
	\describe{
	  \item{\code{ddChangepoint}}{Vector of Date objects corresponding
	    to the changepoints}
	  \item{\code{logLambda}}{Prior on the spline. One of
	    \code{c("iidLogGa","tps","rw1","rw2")}.}
	  \item{\code{tau.gamma}}{}
	  \item{\code{eta.mu}}{}
	  \item{\code{eta.prec}}{}
	  \item{\code{mcmc}}{A names vector of length 3 containing
	    burn-in, number of samples and thinning for the three MCMC chains
	    which are ran. The values are passed on to \code{\link[runjags]{run.jags}}.
	    Default: \code{c(burnin=2500,sample=10000,thin=1)}.}
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
  containing the true number of cases (if possible to compute it based on
  the data in \code{data}. Finally, slot \code{SR} contains the results
  for the proper scoring rules (requires truth to be calculable).
}

\references{
  \enc{Höhle}{Hoehle} M and an der Heiden M (2014), Bayesian Nowcasting
  during the STEC O104:H4 Outbreak in Germany, 2011, Biometrics,
  70(4):993-1002. \doi{10.1111/biom.12194}.\cr
  A preprint is available as
  \url{http://people.su.se/~mhh/pubs/hoehle_anderheiden2014-preprint.pdf}.
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
  \url{http://mcmc-jags.sourceforge.net/}
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
                             mcmc=c(burnin=1000,sample=1000,thin=1))))

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
