\encoding{latin1}
\name{backprojNP}
\alias{backprojNP}
%Internal functions
%\alias{backprojNP.fit}
%\alias{naninf2zero}
%\alias{em.step.becker}

\title{
  Non-parametric back-projection of incidence cases to exposure cases
  using a known incubation time as in Becker et al (1991)
}

\description{
  The function is an implementation of the non-parametric
  back-projection of incidence cases to exposure cases described in
  Becker et al. (1991). The method back-projects exposure times
  from a univariate time series containing the number of symptom onsets per time
  unit. Here, the delay between exposure and symptom onset for an
  individual is seen as a realization of a random variable governed by a
  known probability mass function.
  The back-projection function calculates the expected number of exposures
  \eqn{\lambda_t}{lambda_t} for each time unit under the assumption of a
  Poisson distribution, but without any parametric assumption on how the
  \eqn{\lambda_t}{lambda_t} evolve in time.

  Furthermore, the function contains a bootstrap based procedure, as
  given in Yip et al (2011), which allows an indication of uncertainty
  in the estimated \eqn{\lambda_t}{lambda_T}. The procedure is
  equivalent to the suggestion in Becker and Marschner (1993). However,
  the present implementation in \code{backprojNP} allows only a
  univariate time series, i.e. simultaneous age groups as in Becker and
  Marschner (1993) are not possible.

  The method in Becker et al. (1991) was originally developed for the
  back-projection of AIDS incidence, but it is equally useful for
  analysing the epidemic curve in outbreak situations of a disease
  with long incubation time, e.g. in order to qualitatively investigate
  the effect of intervention measures.
}

\usage{
backprojNP(sts, incu.pmf, 
   control = list(k = 2, 
                  eps = rep(0.005,2), 
                  iter.max=rep(250,2), 
                  Tmark = nrow(sts), 
                  B = -1, 
                  alpha = 0.05, 
                  verbose = FALSE, 
                  lambda0 = NULL,
                  eq3a.method = c("R","C"),
                  hookFun = function(stsbp) {}),
     \dots)
}

\arguments{
  \item{sts}{
    an object of class \code{"\linkS4class{sts}"} (or one that can be
    coerced to that class): contains the observed number of symptom
    onsets as a time series. 
  }
  \item{incu.pmf}{Probability mass function (PMF) of the incubation
    time. The PMF is specified as a vector or matrix with the value of
    the PMF evaluated at \eqn{0,...,d_max}{0,...,d_max}, i.e. note that the
    support includes zero. The value of \eqn{d_max}{d_max} is
    automatically calculated as \code{length(incu.pmf)-1} or
    \code{nrow(incu.pmf)-1}. Note that if the sts object has more
    than one column, then for the backprojection the incubation time is
    either recycled for all components or, if it is a matrix with
    the same number of columns as the sts object, the \eqn{k}{k}'th
    column of \code{incu.pmf} is used for the backprojection of the
    \eqn{k}{k}'th series.
  }
  \item{control}{A list with named arguments controlling the
    functionality of the non-parametric back-projection.
    \describe{
      \item{\code{k}}{An integer representing the smoothing parameter to use in the smoothing step
	of the EMS algorithm. Needs to be an even number.
      }
      \item{\code{eps}}{A vector of length two representing the
	convergence threshold \eqn{\epsilon}{epsilon} of the EMS algorithm, see Details for further
	information. The first value is the threshold to use in the
	\eqn{k=0}{k=0} loop, which forms the values for the parametric
	bootstrap. The second value is the threshold to use in the actual fit
	and bootstrap fitting using the specified \code{k}. If \code{k} is
	only of length one, then this number is replicated twice.
      }
      \item{\code{Tmark}}{Numeric with \eqn{T'\leq T}. Upper time limit
  on which to base convergence, i.e. only the values
  \eqn{\lambda_1,\ldots,\lambda_{T'}} are monitored for convergence. See
  details. 
      }
      \item{\code{iter.max}}{
	The maximum number of EM iterations to do before stopping. 
      }
      \item{\code{B}}{
	Number of parametric bootstrap samples to perform from an
	initial k=0 fit. For each sample a back projection is performed.
	See Becker and Marschner (1993) for details.
      }
      \item{\code{alpha}}{(1-\eqn{\alpha}{alpha})*100\% confidence
  intervals are computed based on the percentile method.
      }
      \item{\code{verbose}}{(boolean). If true show extra progress and
      debug information.
      }
      \item{\code{lambda0}}{Start values for lambda. Vector needs to be
	of the length \code{nrow(sts)}.
      }
      \item{\code{eq3a.method}}{A single character being either
	\code{"R"} or \code{"C"} depending on whether the three nested loops
	of equation 3a in Becker et al. (1991) are to be executed as safe R
	code (can be extremely slow, however the implementation is not
	optimized for speed) or a C code (can be more than 200 times
	faster!). However, the C implementation is experimental and can
	hang R if, e.g., the time series does not go far enough back.
      }
      \item{\code{hookFun}}{
	Hook function called for each iteration of the EM algorithm. The
	function should take a single argument \code{stsbp} of class
	\code{"\linkS4class{stsBP}"} class. It will be have the
        lambda set to the current value of lambda. If no action desired
        just leave the function body empty (default). Additional
	arguments are possible.
      }
    }
  }
  \item{\dots}{Additional arguments are sent to the hook function.
  }
}

\details{
  Becker et al. (1991) specify a non-parametric back-projection
  algorithm based on the Expectation-Maximization-Smoothing (EMS)
  algorithm.
  
  In the present implementation the algorithm iterates until
  \deqn{\frac{||\lambda^{(k+1)} - \lambda^{(k)}||}{||\lambda^{(k)}||} <
  \epsilon} This is a slight adaptation of the proposals in Becker et
  al. (1991). If \eqn{T} is the length of \eqn{\lambda} then one can
  avoid instability of the algorithm near the end by considering only
  the \eqn{\lambda}{lambda}'s with index \eqn{1,\ldots,T'}.

  
  See the references for further information.
}

\value{
  \code{backprojNP} returns an object of \code{"\linkS4class{stsBP}"}.
}

\references{
  Becker NG, Watson LF and Carlin JB (1991), A method for
  non-parametric back-projection and its application to AIDS data,
  Statistics in Medicine, 10:1527-1542.
  
  Becker NG and Marschner IC (1993), A method for estimating the
  age-specific relative risk of HIV infection from AIDS incidence data,
  Biometrika, 80(1):165-178.

  Yip PSF, Lam KF, Xu Y, Chau PH, Xu J, Chang W, Peng Y, Liu Z, Xie X and
  Lau HY (2011), Reconstruction of the Infection Curve for SARS Epidemic
  in Beijing, China Using a Back-Projection Method, Communications in
  Statistics - Simulation and Computation, 37(2):425-433.

  Associations of Age and Sex on Clinical Outcome and Incubation Period
  of Shiga toxin-producing Escherichia coli O104:H4 Infections, 2011
  (2013), Werber D, King LA, \enc{Müller}{Mueller} L, Follin P, Buchholz U, Bernard H,
  Rosner BM, Ethelberg S, de Valk H, \enc{Höhle}{Hoehle} M, American Journal of
  Epidemiology, 178(6):984-992.
}

\author{
  Michael \enc{Höhle}{Hoehle} with help by
  Daniel \enc{Sabanés Bové}{Sabanes Bove} for the \pkg{Rcpp} interface
}

\note{
  The method is still experimental. A proper plot routine for
  \code{stsBP} objects is currently missing.
}

\examples{
#Generate an artificial outbreak of size n starting at time t0 and being of length
n <- 1e3 ; t0 <- 23 ; l <- 10

#PMF of the incubation time is an interval censored gamma distribution
#with mean 15 truncated at 25.
dmax <- 25
inc.pmf <- c(0,(pgamma(1:dmax,15,1.4) - pgamma(0:(dmax-1),15,1.4))/pgamma(dmax,15,1.4))
#Function to sample from the incubation time
rincu <- function(n) {
  sample(0:dmax, size=n, replace=TRUE, prob=inc.pmf)
}
#Sample time of exposure and length of incubation time
set.seed(123)
exposureTimes <- t0 + sample(x=0:(l-1),size=n,replace=TRUE)
symptomTimes <- exposureTimes + rincu(n)

#Time series of exposure (truth) and symptom onset (observed)
X <- table( factor(exposureTimes,levels=1:(max(symptomTimes)+dmax)))
Y <- table( factor(symptomTimes,levels=1:(max(symptomTimes)+dmax)))
#Convert Y to an sts object
Ysts <- sts(Y)

#Plot the outbreak
plot(Ysts, xaxis.labelFormat=NULL, legend=NULL)
#Add true number of exposures to the plot
lines(1:length(Y)+0.2,X,col="red",type="h",lty=2)


#Helper function to show the EM step
plotIt <- function(cur.sts) {
  plot(cur.sts,xaxis.labelFormat=NULL, legend=NULL,ylim=c(0,140))
}

#Call non-parametric back-projection function with hook function but
#without bootstrapped confidence intervals
bpnp.control <- list(k=0,eps=rep(0.005,2),iter.max=rep(250,2),B=-1,hookFun=plotIt,verbose=TRUE)

#Fast C version (use argument: eq3a.method="C")! 
sts.bp <- backprojNP(Ysts, incu.pmf=inc.pmf,
    control=modifyList(bpnp.control,list(eq3a.method="C")), ylim=c(0,max(X,Y)))

#Show result
plot(sts.bp,xaxis.labelFormat=NULL,legend=NULL,lwd=c(1,1,2),lty=c(1,1,1),main="")
lines(1:length(Y)+0.2,X,col="red",type="h",lty=2)

#Do the convolution for the expectation
mu <- matrix(0,ncol=ncol(sts.bp),nrow=nrow(sts.bp))
#Loop over all series
for (j in 1:ncol(sts.bp)) { 
  #Loop over all time points
  for (t in 1:nrow(sts.bp)) {
    #Convolution, note support of inc.pmf starts at zero (move idx by 1)
    i <- seq_len(t)
    mu[t,j] <- sum(inc.pmf[t-i+1] * upperbound(sts.bp)[i,j],na.rm=TRUE)
  }
}
#Show the fit
lines(1:nrow(sts.bp)-0.5,mu[,1],col="green",type="s",lwd=3)

#Non-parametric back-projection including boostrap CIs. B=10 is only
#used for illustration in the documentation example
#In practice use a realistic value of B=1000 or more.
bpnp.control2 <- modifyList(bpnp.control, list(hookFun=NULL,k=2,B=10,eq3a.method="C"))
\dontrun{
bpnp.control2 <- modifyList(bpnp.control, list(hookFun=NULL,k=2,B=1000,eq3a.method="C"))
}
sts.bp2 <- backprojNP(Ysts, incu.pmf=inc.pmf, control=bpnp.control2)

######################################################################
# Plot the result. This is currently a manual routine.
# ToDo: Need to specify a plot method for stsBP objects which also
#       shows the CI.
#
# Parameters:
#  stsBP - object of class stsBP which is to be plotted.
######################################################################

plot.stsBP <- function(stsBP) {
  maxy <- max(observed(stsBP),upperbound(stsBP),stsBP@ci,na.rm=TRUE)
  plot(upperbound(stsBP),type="n",ylim=c(0,maxy), ylab="Cases",xlab="time")
  if (!all(is.na(stsBP@ci))) {
    polygon( c(1:nrow(stsBP),rev(1:nrow(stsBP))),
             c(stsBP@ci[2,,1],rev(stsBP@ci[1,,1])),col="lightgray")
  }
  lines(upperbound(stsBP),type="l",lwd=2)
  legend(x="topright",c(expression(lambda[t])),lty=c(1),col=c(1),fill=c(NA),border=c(NA),lwd=c(2))

  invisible()
}

#Plot the result of k=0 and add truth for comparison. No CIs available
plot.stsBP(sts.bp)
lines(1:length(Y),X,col=2,type="h")
#Same for k=2
plot.stsBP(sts.bp2)
lines(1:length(Y),X,col=2,type="h")

}

\keyword{models}
\keyword{optimize}
