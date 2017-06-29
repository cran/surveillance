\name{findK}
\alias{findK}
\title{Find Reference Value}
\description{
  Calculates the reference value \code{k} for a Poisson or binomial CUSUM
  designed to detect a shift from \eqn{\theta_0} to \eqn{\theta_1}
}
\usage{
findK(theta0, theta1, distr = c("poisson", "binomial"),
      roundK = FALSE, digits = 1, ...)
}
\arguments{
  \item{theta0}{ in-control parameter }
  \item{theta1}{ out-of-control parameter }
  \item{distr}{ \code{"poisson"} or \code{"binomial"} }
  \item{digits}{ the reference value \code{k} is rounded to \code{digits} decimal places}
  \item{roundK}{
    For discrete data and rational reference value there is only
    a limited set of possible values that the CUSUM can take (and
    therefore there is also only a limited set of ARLs).
    If \code{roundK=TRUE}, integer multiples of 0.5 are avoided when
    rounding the reference value \code{k}, %
    i.e. the CUSUM can take more values.}
  \item{\dots}{
    further arguments for the distribution function, i.e. number of
    trials \code{n} for the binomial CDF.}
}
\value{
  Returns reference value \code{k}.
}
\keyword{models}
