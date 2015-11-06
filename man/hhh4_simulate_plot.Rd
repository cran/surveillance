\name{hhh4_simulate_plot}
\alias{plot.hhh4sims}
\alias{as.hhh4simslist}
\alias{plot.hhh4simslist}
\alias{plotHHH4sims_size}
\alias{plotHHH4sims_time}
\alias{scores.hhh4sims}
\alias{scores.hhh4simslist}

\title{
Summarize Simulations from \code{"hhh4"} Models
}

\description{
Arrays of simulated counts from \code{\link{simulate.hhh4}} can be
visualized in various levels of aggregation: final size, time series.
Furthermore, proper scoring rules can be calculated based on the
simulated predictive distributions. Be aware, though, that the current
implementation can only compute univariate scores, i.e., it treats the
predictions at the various time points as independent.
}

\usage{
\method{plot}{hhh4sims}(x, ...)

as.hhh4simslist(x, ...)
\method{plot}{hhh4simslist}(x, type = c("size", "time"), ...,
     groups = NULL, par.settings = list())

plotHHH4sims_size(x, horizontal = TRUE, trafo = NULL, observed = TRUE, ...)

plotHHH4sims_time(x, average = mean, individual = length(x) == 1,
    conf.level = if (individual) 0.95 else NULL,
    matplot.args = list(), initial.args = list(), legend = length(x) > 1,
    xlim = NULL, ylim = NULL, add = FALSE, ...)

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
  \item{type}{
    a character string indicating the summary plot to produce.
  }
  \item{\dots}{
    further arguments passed to methods.
  }
  \item{groups}{
    an optional factor to produce stratified plots by groups of units.
    The special setting \code{groups = TRUE} is a convenient shortcut
    for one plot by unit.
  }
  \item{par.settings}{
    a list of graphical parameters for \code{\link{par}}.
    Sensible defaults for \code{mfrow}, \code{mar} and \code{las} will
    be applied unless overridden or \code{!is.list(par.settings)}.
  }

  \item{horizontal}{
    a logical indicating if the boxplots of the final size distributions
    should be horizontal (the default).
  }
  \item{trafo}{
    an optional transformation function from the \pkg{scales} package, e.g.,
    \code{\link[scales]{sqrt_trans}}.
  }
  \item{observed}{
    a logical indicating if a line and axis value for the observed size
    of the epidemic should be added to the plot.
    Alternatively, a list with graphical parameters can be specified to
    modify the default values.
  }

  \item{average}{
    scalar-valued function to apply to the simulated counts at each time point.
  }
  \item{individual}{
    a logical indicating if the individual simulations should be shown as well.
  }
  \item{conf.level}{
    a scalar in (0,1), which determines the level of the pointwise
    quantiles obtained from the simulated counts at each time point.
    A value of \code{NULL} disables the confidence interval.
  }
  \item{matplot.args}{
    a list of graphical parameters for \code{\link{matlines}}.
  }
  \item{initial.args}{
    if a list (of graphical parameters for \code{\link{lines}}),
    a bar for the initial number of cases is added to the plot.
  }
  \item{legend}{
    a logical or a list of parameters for \code{\link{legend}}.
  }
  \item{xlim,ylim}{
    vectors of length 2 determining the axis limits.
  }
  \item{add}{
    a logical indicating if the (mean) simulated time series should be
    added to an existing plot.
  }

  \item{which}{
    a character vector indicating which proper scoring rules to compute.
    By default, only the ranked probability score (\code{"rps"}) is
    calculated. Other options include \code{"logs"} and \code{"dss"}.
  }
  \item{units}{
    if non-\code{NULL}, an integer or character vector indexing the
    columns of \code{x} to compute the scores for that subset only.
  }
  \item{drop}{
    a logical indicating if univariate dimensions should be dropped
    (the default).
  }
}

%% \value{
%% }

\author{
  Sebastian Meyer
}

\examples{
### univariate example
data("salmAllOnset")

## fit a hhh4 model to the first 13 years
salmModel <- list(end = list(f = addSeason2formula(~1 + t)),
                  ar = list(f = ~1), family = "NegBin1", subset = 2:678)
salmFit <- hhh4(salmAllOnset, salmModel)

## simulate the next 20 weeks ahead
salmSims <- simulate(salmFit, nsim = 300, seed = 3, subset = 678 + seq_len(20),
                     y.start = observed(salmAllOnset)[678,])

## compare final size distribution to observed value
plot(salmSims)

## simulated time series
plot(salmSims, type = "time", main = "2-weeks-ahead simulation")


### multivariate example
data("measlesWeserEms")

## fit a hhh4 model to the first year
measlesModel <- list(
    end = list(f = addSeason2formula(~1), offset = population(measlesWeserEms)),
    ar = list(f = ~1),
    ne = list(f = ~1 + log(pop),
        weights = W_powerlaw(maxlag = 5, normalize = TRUE)),
    family = "NegBin1", subset = 2:52,
    data = list(pop = population(measlesWeserEms)))
measlesFit1 <- hhh4(measlesWeserEms, control = measlesModel)
measlesFit2 <- update(measlesFit1, family = "Poisson")

## simulate realizations from this model during the second year
measlesSims <- lapply(X = list(NegBin = measlesFit1, Poisson = measlesFit2),
                      FUN = simulate, nsim = 50, seed = 1, subset = 53:104,
                      y.start = observed(measlesWeserEms)[52,])

## final size of the first model
plot(measlesSims[[1]])

## stratified by groups of districts
plot(measlesSims[[1]], groups = factor(substr(colnames(measlesWeserEms), 4, 4)))

## a class and plot-method for a list of simulations from different models
measlesSims <- as.hhh4simslist(measlesSims)
plot(measlesSims)

## simulated time series
plot(measlesSims, type = "time", individual = TRUE, ylim = c(0, 80))

## compare proper scoring rules for a specific subset of the regions
## (CAVE: these are univariate scores for each time point and region,
##        which do not account for dependence over time)
measlesScores5 <- scores(measlesSims, which = "rps",
                         units = substr(colnames(measlesWeserEms), 4, 4) == "5")
sapply(measlesScores5, mean)
}

\keyword{hplot}
