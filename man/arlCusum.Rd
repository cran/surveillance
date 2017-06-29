\name{arlCusum}
\alias{arlCusum}
\title{Calculation of Average Run Length for discrete CUSUM schemes}
\description{
  Calculates the average run length (ARL) for an upward CUSUM scheme for discrete
  distributions (i.e. Poisson and binomial) using the Markov chain approach.
}
\usage{
arlCusum(h=10, k=3, theta=2.4, distr=c("poisson","binomial"),
         W=NULL, digits=1, ...)
}
\arguments{
  \item{h}{ decision interval}
  \item{k}{ reference value}
  \item{theta}{distribution parameter for the cumulative distribution function
               (cdf) \eqn{F}, i.e. rate
               \eqn{\lambda} for Poisson variates or probability \eqn{p}
               for binomial variates}
  \item{distr}{ \code{"poisson"} or \code{"binomial"}
               }      %ppois, pbinom
  \item{W}{Winsorizing value \code{W} for a robust CUSUM,
          to get a nonrobust CUSUM set %\code{W} is set  to
          \code{W} > \code{k}+\code{h}. If \code{NULL}, a nonrobust CUSUM is used.}
  \item{digits}{ \code{k} and \code{h} are rounded to \code{digits} decimal places }
  \item{\dots}{ further arguments for the distribution function, i.e. number of trials \code{n}
                for binomial cdf }
}
\value{
Returns a list with the ARL of the regular (zero-start)
and the fast initial response (FIR)
CUSUM scheme with reference value \code{k}, decision interval \code{h} for
\eqn{X \sim F(\theta)}, where F is the Poisson or binomial CDF.
\item{ARL}{one-sided ARL of the regular (zero-start) CUSUM scheme}
\item{FIR.ARL}{one-sided ARL of the FIR CUSUM scheme with head start
  \eqn{\frac{\code{h}}{2}} }
}
\keyword{models}
\source{Based on the FORTRAN code of

Hawkins, D. M. (1992). Evaluation of Average Run Lengths of Cumulative
  Sum Charts for an Arbitrary Data Distribution. Communications in
  Statistics - Simulation and Computation, 21(4), p. 1001-1020.
}
