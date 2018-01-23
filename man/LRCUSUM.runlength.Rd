\name{LRCUSUM.runlength}
\alias{LRCUSUM.runlength}
\alias{outcomeFunStandard}
\alias{LLR.fun}
\encoding{latin1}

\title{Run length computation of a CUSUM detector}
\description{
  Compute run length for a count data or categorical CUSUM. The
  computations are based on a Markov representation of the
  likelihood ratio based CUSUM.
}
\usage{
LRCUSUM.runlength(mu,mu0,mu1,h,dfun, n, g=5,outcomeFun=NULL,...)
}
\arguments{
  \item{mu}{\eqn{k-1 \times T} matrix with true proportions, i.e. equal to mu0 or mu1 if one wants to compute e.g. \eqn{ARL_0} or \eqn{ARL_1}.}
  \item{mu0}{\eqn{k-1 \times T} matrix with in-control proportions}
  \item{mu1}{\eqn{k-1 \times T} matrix with out-of-control proportion}
  \item{h}{The threshold h which is used for the CUSUM.}
  \item{dfun}{The probability mass function or density used to compute
    the likelihood ratios of the CUSUM. In a negative binomial CUSUM
    this is \code{dnbinom}, in a binomial CUSUM \code{dbinom} and in a
    multinomial CUSUM \code{dmultinom}.}
  \item{n}{Vector of length \eqn{T} containing the total number of
    experiments for each time point.}
  \item{g}{The number of levels to cut the state space into when
    performing the Markov chain approximation. Sometimes also denoted
    \eqn{M}. Note that the quality of the approximation depends very
    much on \eqn{g}. If \eqn{T} greater than, say, 50 its necessary to
  increase the value of \eqn{g}.}
  \item{outcomeFun}{A hook function to compute all possible outcome
    states to compute the likelihood ratio for. If \code{NULL} then the
    default function \code{outcomeFunStandard(k,n)} is used. This
    function uses the Cartesian product of \code{0:n} for \code{k} components.}
  \item{\dots}{Additional arguments to send to \code{dfun}.}
}
\details{
Brook and Evans (1972) formulated an approximate approach based
on Markov chains to determine the PMF of the run length of a
time-constant CUSUM detector. They describe the dynamics of the CUSUM
statistic by a Markov chain with a discretized state space of
size \eqn{g+2}. This is adopted to the time varying case in
\enc{Höhle}{Hoehle} (2010) and implemented in R using the \dots notation
such that it works for a very large class of distributions.

}
\seealso{\code{\link{categoricalCUSUM}}}
\value{A list with five components
  \item{P}{An array of \eqn{g+2 \times g+2} transition matrices of the
    approximation Markov chain.}
  \item{pmf}{Probability mass function (up to length \eqn{T}) of the run
    length variable.}
  \item{cdf}{Cumulative density function (up to length \eqn{T}) of the run
    length variable.}
  \item{arl}{If the model is time homogenous (i.e. if \eqn{T==1}) then
    the ARL is computed based on the stationary distribution of the
    Markov chain. See the eqns in the reference for details. Note: If
    the model is not time homogeneous then the function returns
    \code{NA} and the ARL has to be approximated manually from the
    output. One could use \code{sum(1:length(pmf) * pmf)}, which is an
    approximation because of using a finite support for a sum which
    should be from 1 to infinity.
  }
}
\references{
  \enc{Höhle}{Hoehle}, M. (2010):
  Online change-point detection in categorical time series. 
  In: T. Kneib and G. Tutz (Eds.), Statistical
  Modelling and Regression Structures - Festschrift in Honour of Ludwig
  Fahrmeir, Physica-Verlag, pp. 377-397. Preprint available as
  \url{http://www.math.su.se/~hoehle/pubs/hoehle2010-preprint.pdf}

  \enc{Höhle}{Hoehle}, M. and Mazick, A. (2010):
  Aberration detection in R illustrated
  by Danish mortality monitoring. In: T. Kass-Hout
  and X. Zhang (Eds.), Biosurveillance: A Health Protection Priority,
  CRCPress. Preprint available as
  \url{http://www.math.su.se/~hoehle/pubs/hoehle_mazick2009-preprint.pdf}
  
  Brook, D. and Evans, D. A. (1972),  An approach to the probability
  distribution of Cusum run length, Biometrika,
  59:3, pp. 539--549.
}
\examples{
######################################################
#Run length of a time constant negative binomial CUSUM
######################################################

#In-control and out of control parameters
mu0 <- 10
alpha <- 1/2
kappa <- 2

#Density for comparison in the negative binomial distribution
dY <- function(y,mu,log=FALSE, alpha, ...) {
  dnbinom(y, mu=mu, size=1/alpha, log=log)
}

#In this case "n" is the maximum value to investigate the LLR for
#It is assumed that beyond n the LLR is too unlikely to be worth
#computing.
LRCUSUM.runlength( mu=t(mu0), mu0=t(mu0), mu1=kappa*t(mu0), h=5,
  dfun = dY, n=rep(100,length(mu0)), alpha=alpha)

h.grid <- seq(3,6,by=0.3)
arls <- sapply(h.grid, function(h) {
  LRCUSUM.runlength( mu=t(mu0), mu0=t(mu0), mu1=kappa*t(mu0), h=h,
  dfun = dY, n=rep(100,length(mu0)), alpha=alpha,g=20)$arl
})
plot(h.grid, arls,type="l",xlab="threshold h",ylab=expression(ARL[0]))


if (surveillance.options("allExamples"))
{
  ######################################################
  #Run length of a time varying negative binomial CUSUM
  ######################################################

  mu0 <- matrix(5*sin(2*pi/52 * 1:200) + 10,ncol=1)
  
  rl <- LRCUSUM.runlength( mu=t(mu0), mu0=t(mu0), mu1=kappa*t(mu0), h=2,
    dfun = dY, n=rep(100,length(mu0)), alpha=alpha,g=20)

  plot(1:length(mu0),rl$pmf,type="l",xlab="t",ylab="PMF")
  plot(1:length(mu0),rl$cdf,type="l",xlab="t",ylab="CDF")
}


########################################################
# Further examples contain the binomial, beta-binomial
# and multinomial CUSUMs. Hopefully, these will be added
# in the future.
########################################################

#dfun function for the multinomial distribution (Note: Only k-1 categories are specified).
dmult <- function(y, size,mu, log = FALSE) {
    return(dmultinom(c(y,size-sum(y)), size = size, prob=c(mu,1-sum(mu)), log = log))
}

#Example for the time-constant multinomial distribution
#with size 100 and in-control and out-of-control parameters as below.
n <- 100
pi0 <- as.matrix(c(0.5,0.3,0.2))
pi1 <- as.matrix(c(0.38,0.46,0.16))

#ARL_0
LRCUSUM.runlength(mu=pi0[1:2,,drop=FALSE],mu0=pi0[1:2,,drop=FALSE],mu1=pi1[1:2,,drop=FALSE],
                  h=5,dfun=dmult, n=n, g=15)$arl
#ARL_1
LRCUSUM.runlength(mu=pi1[1:2,,drop=FALSE],mu0=pi0[1:2,,drop=FALSE],mu1=pi1[1:2,,drop=FALSE],
                  h=5,dfun=dmult, n=n, g=15)$arl


}

\author{M. \enc{Höhle}{Hoehle}}
\keyword{regression}

