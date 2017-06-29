\name{findH}
\alias{findH}
\alias{hValues}
\title{Find decision interval for given in-control ARL and reference value}
\description{
  Function to find a decision interval \code{h}* for given reference value \code{k}
  and desired ARL \eqn{\gamma} so that the
  average run length for a Poisson or Binomial CUSUM with in-control
  parameter \eqn{\theta_0}, reference value \code{k} and is approximately \eqn{\gamma},
  i.e. \eqn{\Big| \frac{ARL(h^*) -\gamma}{\gamma} \Big| < \epsilon},
  or larger, i.e.
  \eqn{ARL(h^*) > \gamma }.
}
\usage{
findH(ARL0, theta0, s = 1, rel.tol = 0.03, roundK = TRUE,
      distr = c("poisson", "binomial"), digits = 1, FIR = FALSE, ...)

hValues(theta0, ARL0, rel.tol=0.02, s = 1, roundK = TRUE, digits = 1,
        distr = c("poisson", "binomial"), FIR = FALSE, ...)
}
\arguments{
  \item{ARL0}{ desired in-control ARL \eqn{\gamma} }
  \item{theta0}{in-control parameter \eqn{\theta_0}}
  \item{s}{change to detect, see details}
  \item{distr}{ \code{"poisson"} or \code{"binomial"} }
  \item{rel.tol}{relative tolerance, i.e. the search for \code{h}* is
                 stopped if \eqn{\Big| \frac{ARL(h^*) -\gamma}{\gamma} \Big| < } \code{rel.tol} }
  \item{digits}{the reference value \code{k} and the decision  interval \code{h}
                are rounded to \code{digits} decimal places}
  \item{roundK}{ passed to \code{findK} }
  \item{FIR}{if \code{TRUE}, the decision interval that leads to the desired ARL
                 for a FIR CUSUM with head start
                 \eqn{\frac{\code{h}}{2}} is returned }
  \item{\dots}{ further arguments for the distribution function, i.e. number
                of trials \code{n} for binomial cdf }
}
\value{
  \code{findH} returns a vector and \code{hValues} returns a matrix with elements
  \item{theta0}{in-control parameter}
  \item{h}{decision interval}
  \item{k}{reference value}
  \item{ARL}{ARL for a CUSUM with parameters \code{k} and \code{h} }
  \item{rel.tol}{corresponds to \eqn{\Big| \frac{ARL(h) -\gamma}{\gamma} \Big|} }
}
\details{
The out-of-control parameter used to determine the reference value \code{k}
is specified as:
\deqn{\theta_1 = \lambda_0 + s \sqrt{\lambda_0} }
for a Poisson variate \eqn{X \sim Po(\lambda)}

\deqn{\theta_1 = \frac{s \pi_0}{1+(s-1) \pi_0} }
for a Binomial variate \eqn{X \sim Bin(n, \pi) }
}
\keyword{models}
