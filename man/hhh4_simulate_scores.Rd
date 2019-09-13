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

## simulate the next 20 weeks ahead (with very small 'nsim' for speed)
salmSims <- simulate(salmFit, nsim = 500, seed = 3, subset = 678 + seq_len(20),
                     y.start = observed(salmAllOnset)[678,])
if (requireNamespace("fanplot"))
    plot(salmSims, "fan")


### calculate scores at each time point

## using empirical distribution of simulated counts as forecast distribution
scores(salmSims, which = c("rps", "logs", "dss"))
## observed count sometimes not covered by simulations -> infinite log-score
## => for a more detailed forecast, either considerably increase 'nsim', or:

## 1. use continuous density() of simulated counts as forecast distribution
fi <- apply(salmSims, 1, function (x) approxfun(density(x)))
logs_kde <- mapply(function (f, y) -log(f(y)),
                   f = fi, y = observed(attr(salmSims,"stsObserved")))
cbind("empirical" = scores(salmSims, "logs"), "density" = logs_kde)
## a similar KDE approach is implemented in scoringRules::logs_sample()

## 2. average conditional predictive NegBin's of simulated trajectories,
##    currently only implemented in HIDDA.forecasting::dhhh4sims()


\dontrun{
### produce a PIT histogram

## using empirical distribution of simulated counts as forecast distribition
pit(x = observed(attr(salmSims, "stsObserved")),
    pdistr = apply(salmSims, 1:2, ecdf))
## long-term forecast is badly calibrated (lower tail is unused, see fan above)
## we also get a warning for the same reason as infinite log-scores
}
}

\keyword{univar}
