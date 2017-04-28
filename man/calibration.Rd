\name{calibrationTest}
\alias{calibrationTest}
\alias{calibrationTest.default}

\title{
  Calibration Tests for Poisson or Negative Binomial Predictions
}

\description{
  The implemented calibration tests for Poisson or negative binomial
  predictions of count data are based on proper scoring rules and
  described in detail in Wei and Held (2014).
  The following proper scoring rules are available:
  Dawid-Sebastiani score (\code{"dss"}),
  logarithmic score (\code{"logs"}),
  ranked probability score (\code{"rps"}).
}

\usage{
calibrationTest(x, ...)

\method{calibrationTest}{default}(x, mu, size = NULL,
                which = c("dss", "logs", "rps"),
                tolerance = 1e-4, method = 2, ...)
}

\arguments{
  \item{x}{
    the observed counts.
    All involved functions are vectorized and also accept matrices or arrays.
  }
  \item{mu}{
    the means of the predictive distributions for the
    observations \code{x}.
  }
  \item{size}{
    either \code{NULL} (default), indicating Poisson predictions with mean
    \code{mu}, or dispersion parameters of
    negative binomial forecasts for the observations \code{x},
    parametrized as in \code{\link{dnbinom}} with variance
    \code{mu*(1+mu/size)}.
  }
  \item{which}{
    a character string indicating which proper scoring rule to apply.
  }
  \item{tolerance}{
    absolute tolerance for the null expectation and variance of
    \code{"logs"} and \code{"rps"}. For the latter, see the note below.
    Unused for \code{which = "dss"} (closed form).
  }
  \item{method}{
    selection of the \eqn{z}-statistic: \code{method = 2} refers to
    the alternative test statistic \eqn{Z_s^*} of Wei and Held (2014,
    Discussion), which has been recommended for low counts.
    \code{method = 1} corresponds to Equation 5 in Wei and Held (2014).
  }
  \item{\dots}{
    unused (argument of the generic).
  }
}

\value{
  an object of class \code{"htest"},
  which is a list with the following components:
  \item{method}{a character string indicating the type of test
    performed (including \code{which} scoring rule).}
  \item{data.name}{a character string naming the supplied \code{x} argument.}
  \item{statistic}{the \eqn{z}-statistic of the test.}
  \item{parameter}{the number of predictions underlying the test, i.e., \code{length(x)}.}
  \item{p.value}{the p-value for the test.}
}

\note{
  If the \CRANpkg{gsl} package is installed, its implementations of the
  Bessel and hypergeometric functions are used when calculating the null
  expectation and variance of the \code{rps}.
  These functions are faster and yield more accurate results (especially
  for larger \code{mu}).
}

\references{
  Wei, W. and Held, L. (2014):
  Calibration tests for count data.
  \emph{Test}, \bold{23}, 787-805.  
}

\author{
  Sebastian Meyer and Wei Wei
}

\examples{
mu <- c(0.1, 1, 3, 6, pi, 100)
size <- 0.1
set.seed(1)
y <- rnbinom(length(mu), mu = mu, size = size)
calibrationTest(y, mu = mu, size = size) # p = 0.99
calibrationTest(y, mu = mu, size = 1) # p = 4.3e-05
calibrationTest(y, mu = 1, size = size) # p = 0.6959
calibrationTest(y, mu = 1, size = size, which = "rps") # p = 0.1286
}

\keyword{htest}
