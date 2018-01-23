\name{hhh4_simulate_scores}
\alias{scores.hhh4sims}
\alias{scores.hhh4simslist}

\title{
  Proper Scoring Rules for Simulations from \code{hhh4} Models
}

\description{
Calculate proper scoring rules based on simulated predictive distributions.
}

\usage{
\method{scores}{hhh4sims}(x, which = "rps", units = NULL, ..., drop = TRUE)
\method{scores}{hhh4simslist}(x, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"hhh4sims"} (as resulting from the
    \code{\link[=simulate.hhh4]{simulate}}-method for
    \code{"\link{hhh4}"} models if \code{simplify = TRUE} was set),
    or an \code{"hhh4simslist"}, i.e.,
    a list of such simulations potentially obtained from different
    model fits (using the same simulation period).
  }
  \item{which}{
    a character vector indicating which proper scoring rules to compute.
    By default, only the ranked probability score (\code{"rps"}) is
    calculated. Other options include \code{"logs"} and \code{"dss"}.
  }
  \item{units}{
    if non-\code{NULL}, an integer or character vector indexing the
    columns of \code{x} for which to compute the scores.
  }
  \item{drop}{
    a logical indicating if univariate dimensions should be dropped
    (the default).
  }
  \item{\dots}{
    unused (argument of the generic).
  }
}

\details{
  This implementation can only compute \emph{univariate scores}, i.e.,
  independently for each time point.
  
  The logarithmic score is badly estimated if the domain is large and
  there are not enough samples to cover the underlying distribution in
  enough detail (the score becomes infinite when an observed value does
  not occur in the samples). An alternative is to use kernel density
  estimation as implemented in \code{\link[scoringRules]{logs_sample}}
  in the \R package \CRANpkg{scoringRules}.
}

\author{
  Sebastian Meyer
}

\examples{
data("salmAllOnset")

## fit a hhh4 model to the first 13 years
salmModel <- list(end = list(f = addSeason2formula(~1 + t)),
                  ar = list(f = ~1), family = "NegBin1", subset = 2:678)
salmFit <- hhh4(salmAllOnset, salmModel)

## simulate the next 20 weeks ahead
salmSims <- simulate(salmFit, nsim = 500, seed = 3, subset = 678 + seq_len(20),
                     y.start = observed(salmAllOnset)[678,])

## calculate the RPS at each time point
scores(salmSims, which = "rps")

## produce a PIT histogram based on the empirical distribution function
## of the simulated counts as the forecast distribution at each time point
pit(x = observed(attr(salmSims, "stsObserved")),
    pdistr = apply(salmSims, 1:2, ecdf))
}

\keyword{univar}
