\name{siaf.simulatePC}
\alias{siaf.simulatePC}
\title{
  Simulation from an Isotropic Spatial Kernel via Polar Coordinates
}
\description{
  To sample points from isotropic spatial kernels
  \eqn{f_2(s) = f(||s||)} such as \code{\link{siaf.powerlaw}} on a
  bounded domain (i.e., \eqn{||s|| < \code{ub}}), it is
  convenient to switch to polar coordinates \eqn{(r,\theta)},
  which have a density proportional to
  \eqn{r f_2((r \cos(\theta), r \sin(\theta))) = r f(r)}
  (independent of the angle \eqn{\theta} due to isotropy).
  The angle is thus simply drawn uniformly in \eqn{[0,2\pi)}, and
  \eqn{r} can be sampled by the inversion method, where numeric root
  finding is used for the quantiles (since the quantile function is not
  available in closed form).
}

\usage{
siaf.simulatePC(intrfr)
}

\arguments{
  \item{intrfr}{
    a function computing the integral of \eqn{r f(r)} from 0 to \code{R}
    (first argument, not necessarily named \code{R}). Parameters of the
    function are passed as its second argument and a third argument is
    the event type.
  }
}

\value{
  a function with arguments \code{(n, siafpars, type, ub)}, which
  samples \code{n} points from the spatial kernel \eqn{f_2(s)} within the
  disc of radius \code{ub}, where \code{siafpars} and \code{type} are
  passed as second and third argument to \code{intrfr}.
  The environment of the returned function will be the caller's environment.
}

\author{
  Sebastian Meyer
}

\examples{
simfun <- siaf.powerlaw()$simulate
## is internally generated as siaf.simulatePC(intrfr.powerlaw)

set.seed(1)
simfun(n=10, siafpars=log(c(sigma=1, d=2)), ub=5)
}

\keyword{internal}
