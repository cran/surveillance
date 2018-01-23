\name{hhh4_simulate_plot}
\alias{plot.hhh4sims}
\alias{aggregate.hhh4sims}
\alias{as.hhh4simslist}
\alias{plot.hhh4simslist}
\alias{aggregate.hhh4simslist}
\alias{plotHHH4sims_size}
\alias{plotHHH4sims_time}
\alias{plotHHH4sims_fan}

\title{
Plot Simulations from \code{"hhh4"} Models
}

\description{
Arrays of simulated counts from \code{\link{simulate.hhh4}} can be
visualized as final size boxplots, individual or average time series,
or fan charts (using the \CRANpkg{fanplot} package).
An \code{aggregate}-method is also available.
}

\usage{
\method{plot}{hhh4sims}(x, ...)
\method{aggregate}{hhh4sims}(x, units = TRUE, time = FALSE, ..., drop = FALSE)

as.hhh4simslist(x, ...)
\method{plot}{hhh4simslist}(x, type = c("size", "time", "fan"), ...,
     groups = NULL, par.settings = list())
\method{aggregate}{hhh4simslist}(x, units = TRUE, time = FALSE, ..., drop = FALSE)

plotHHH4sims_size(x, horizontal = TRUE, trafo = NULL, observed = TRUE,
                  names = base::names(x), ...)

plotHHH4sims_time(x, average = mean, individual = length(x) == 1,
    conf.level = if (individual) 0.95 else NULL,
    matplot.args = list(), initial.args = list(), legend = length(x) > 1,
    xlim = NULL, ylim = NULL, add = FALSE, ...)

plotHHH4sims_fan(x, which = 1,
    fan.args = list(), observed.args = list(), initial.args = list(),
    means.args = NULL, key.args = NULL, xlim = NULL, ylim = NULL,
    add = FALSE, xaxis = list(), ...)
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
  \item{names}{
    a character vector of names for \code{x}.
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
    a logical, a character vector (providing names for \code{x}),
    or a list of parameters for \code{\link{legend}}.
  }
  \item{xlim,ylim}{
    vectors of length 2 determining the axis limits.
  }
  \item{add}{
    a logical indicating if the (mean) simulated time series or the fan
    chart, respectively, should be added to an existing plot.
  }
  \item{which}{
    a single integer or a character string selecting the model in
    \code{x} for which to produce the fan chart.
    This is only relevant if \code{x} is a \code{"hhh4simslist"} of
    simulations from multiple models.
    Defaults to the first model.
  }
  \item{fan.args}{
    a list of graphical parameters for the \code{\link[fanplot]{fan}},
    e.g., to employ a different \code{\link{colorRampPalette}} as
    \code{fan.col}, or to enable contour lines via \code{ln}.
  }
  \item{observed.args}{
    if a list (of graphical parameters for \code{\link{lines}}),
    the originally observed counts are added to the plot.
  }
  \item{means.args}{
    if a list (of graphical parameters for \code{\link{lines}}),
    the point forecasts are added to the plot (by default as
    a white line within the fan).
  }
  \item{key.args}{
    if a list, a color key (in \code{\link[fanplot]{fan}}'s
    \code{"boxfan"}-style) is added to the fan chart. The list may
    include positioning parameters \code{start} (the x-position) and
    \code{ylim} (the y-range of the color key), \code{space} to modify
    the width of the boxfan, and \code{rlab} to modify the labels.
    The color key is disabled by default.
    An alternative way of labeling the quantiles is via the argument
    \code{ln} in \code{fan.args}, see the Examples.
  }
  \item{xaxis}{
    if a list of arguments for \code{\link{addFormattedXAxis}}, that
    function is used to draw the time axis, otherwise a default x-axis
    is drawn.
  }

  \item{units}{
    a logical indicating aggregation over units. Can also be a factor
    (or something convertible to a factor using \code{\link{as.factor}})
    to aggregate groups of units.
  }
  \item{time}{
    a logical indicating if the counts should be summed over the whole
    simulation period.
  }
  \item{drop}{
    a logical indicating if the unit dimension and the \code{"hhh4sims"}
    (or \code{"hhh4simslist"}) class should be dropped after aggregating
    over (groups of) units.
  }
}

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
summary(aggregate(salmSims, time = TRUE))  # summary of simulated values
plot(salmSims, type = "size")

## individual and average simulated time series with a confidence interval
plot(salmSims, type = "time", main = "20-weeks-ahead simulation")

## fan chart based on the quantiles of the simulated counts at each time point
## point forecasts are represented by a white line within the fan
if (requireNamespace("fanplot")) {
    plot(salmSims, type = "fan", main = "20-weeks-ahead simulation",
         fan.args = list(ln = 1:9/10), means.args = list())
}


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

## use a Poisson distribution instead (just for comparison)
measlesFit2 <- update(measlesFit1, family = "Poisson")

## simulate realizations from these models during the second year
measlesSims <- lapply(X = list(NegBin = measlesFit1, Poisson = measlesFit2),
                      FUN = simulate, nsim = 50, seed = 1, subset = 53:104,
                      y.start = observed(measlesWeserEms)[52,])

## final size of the first model
plot(measlesSims[[1]])

## stratified by groups of districts
mygroups <- factor(substr(colnames(measlesWeserEms), 4, 4))
apply(aggregate(measlesSims[[1]], time = TRUE, units = mygroups), 1, summary)
plot(measlesSims[[1]], groups = mygroups)

## a class and plot-method for a list of simulations from different models
measlesSims <- as.hhh4simslist(measlesSims)
plot(measlesSims)

## simulated time series
plot(measlesSims, type = "time", individual = TRUE, ylim = c(0, 80))

## fan charts
if (requireNamespace("fanplot")) {
    opar <- par(mfrow = c(2,1))
    plot(measlesSims, type = "fan", which = 1, ylim = c(0, 80), main = "NegBin",
         key.args = list())
    plot(measlesSims, type = "fan", which = 2, ylim = c(0, 80), main = "Poisson")
    par(opar)
}
}

\keyword{hplot}
