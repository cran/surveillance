\name{scores}
\alias{scores}
\alias{scores.default}
\alias{logs}
\alias{rps}
\alias{dss}
\alias{ses}

\title{
  Proper Scoring Rules for Poisson or Negative Binomial Predictions
}

\description{
  Proper scoring rules for Poisson or negative binomial predictions
  of count data are described in Czado et al. (2009).
  The following scores are implemented:
  logarithmic score (\code{logs}),
  ranked probability score (\code{rps}),
  Dawid-Sebastiani score (\code{dss}),
  squared error score (\code{ses}).
}

\usage{
scores(x, ...)

\method{scores}{default}(x, mu, size = NULL,
       which = c("logs", "rps", "dss", "ses"),
       sign = FALSE, ...)

logs(x, mu, size = NULL)
rps(x, mu, size = NULL, k = 40, tolerance = sqrt(.Machine$double.eps))
dss(x, mu, size = NULL)
ses(x, mu, size = NULL)
}

\arguments{
  \item{x}{
    the observed counts.
    All functions are vectorized and also accept matrices or arrays.
    Dimensions are preserved.
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
    a character vector specifying which scoring rules to apply.
    By default, all four proper scores are calculated.
    The normalized squared error score (\code{"nses"}) is also available
    but it is improper and hence not computed by default.
  }
  \item{sign}{
    a logical indicating if the function should also return
    \code{sign(x-mu)}, i.e., the sign of the difference between
    the observed counts and corresponding predictions.
  }
  \item{\dots}{
    unused (argument of the generic).
  }
  \item{k}{
    scalar argument controlling the finite sum approximation for the
    \code{rps} with truncation at \code{ceiling(mu + k*sd)}.
  }
  \item{tolerance}{
    absolute tolerance for the finite sum approximation employed in the
    \code{rps} calculation. A warning is produced if the approximation
    with \code{k} summands is insufficient for the specified
    \code{tolerance}. In this case, increase \code{k} for higher
    precision (or use a larger tolerance).
  }
}

\value{
  The scoring functions return the individual scores for the predictions
  of the observations in \code{x} (maintaining their dimension attributes).
  
  The default \code{scores}-method applies the selected (\code{which})
  scoring functions (and calculates \code{sign(x-mu)}) and returns the
  results in an array (via \code{\link{simplify2array}}), where the last
  dimension corresponds to the different scores.
}

\references{
  Czado, C., Gneiting, T. and Held, L. (2009):
    Predictive model assessment for count data.
    \emph{Biometrics}, \bold{65} (4), 1254-1261.
    \doi{10.1111/j.1541-0420.2009.01191.x}
}

\seealso{
  The R package \CRANpkg{scoringRules} implements the logarithmic
  score and the (continuous) ranked probability score for many
  distributions.
}

\author{
  Sebastian Meyer and Michaela Paul
}

\examples{
mu <- c(0.1, 1, 3, 6, pi, 100)
size <- 0.1
set.seed(1)
y <- rnbinom(length(mu), mu = mu, size = size)
scores(y, mu = mu, size = size)
scores(y, mu = mu, size = 1)  # ses ignores the variance
scores(y, mu = 1, size = size)

## apply a specific scoring rule
scores(y, mu = mu, size = size, which = "rps")
rps(y, mu = mu, size = size)
}

\keyword{univar}
