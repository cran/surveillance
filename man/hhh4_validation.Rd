\name{hhh4_validation}
\alias{oneStepAhead}
\alias{quantile.oneStepAhead}
\alias{confint.oneStepAhead}
\alias{plot.oneStepAhead}
\alias{scores.oneStepAhead}
\alias{scores.hhh4}
\alias{calibrationTest.oneStepAhead}
\alias{calibrationTest.hhh4}
\alias{pit.oneStepAhead}
\alias{pit.hhh4}

\title{Predictive Model Assessment for \code{hhh4} Models}
\description{
  The function \code{oneStepAhead} computes successive one-step-ahead
  predictions for a (random effects) HHH model fitted by \code{\link{hhh4}}.
  These can be inspected using the \code{quantile}, \code{confint} or
  \code{plot} methods.
  The associated \code{\link{scores}}-method computes a number of (strictly) proper
  scoring rules based on such one-step-ahead predictions;
  see Paul and Held (2011) for details.
  There are also \code{\link{calibrationTest}} and \code{\link{pit}}
  methods for \code{oneStepAhead} predictions.
  
  Scores, calibration tests and PIT histograms can also be
  computed for the fitted values of an \code{hhh4} model
  (i.e., in-sample/training data evaluation).
}

\usage{
oneStepAhead(result, tp, type = c("rolling", "first", "final"),
             which.start = c("current", "final"),
             keep.estimates = FALSE, verbose = TRUE, cores = 1)

\method{quantile}{oneStepAhead}(x, probs = c(2.5, 10, 50, 90, 97.5)/100, ...)
\method{confint}{oneStepAhead}(object, parm, level = 0.95, ...)
\method{plot}{oneStepAhead}(x, unit = 1, probs = 1:99/100,
     start = NULL, means.args = NULL, ...)

## assessment of "oneStepAhead" predictions
\method{scores}{oneStepAhead}(x, which = c("logs", "rps", "dss", "ses"),
       units = NULL, sign = FALSE, individual = FALSE, reverse = FALSE, ...)
\method{calibrationTest}{oneStepAhead}(x, units = NULL, ...)
\method{pit}{oneStepAhead}(x, units = NULL, ...)

## assessment of the "hhh4" model fit (in-sample predictions)
\method{scores}{hhh4}(x, which = c("logs", "rps", "dss", "ses"),
       subset = x$control$subset, units = seq_len(x$nUnit), sign = FALSE, ...)
\method{calibrationTest}{hhh4}(x,
                subset = x$control$subset, units = seq_len(x$nUnit), ...)
\method{pit}{hhh4}(x, subset = x$control$subset, units = seq_len(x$nUnit), ...)
}

\arguments{
  \item{result}{fitted \code{\link{hhh4}} model (class \code{"hhh4"}).}
  \item{tp}{
    numeric vector of length 2 specifying the time range in
    which to compute one-step-ahead predictions (for the time points
    \code{tp[1]+1}, \ldots, \code{tp[2]+1}).
    If a single time index is specified, it is interpreted as
    \code{tp[1]}, and \code{tp[2]} is set to the penultimate time point
    of \code{result$control$subset}.
  }
  \item{type}{
    The default \code{"rolling"} procedure sequentially 
    refits the model up to each time point in \code{tp} and computes
    the one-step-ahead predictions for the respective next time point.
    The alternative \code{type}s are no true one-step-ahead predictions
    but much faster:
    \code{"first"} will refit the model for the first time point
    \code{tp[1]} only and use this specific fit to calculate all
    subsequent predictions, whereas
    \code{"final"} will just use \code{result} to calculate these.
    The latter case thus gives nothing else than a subset of
    \code{result$fitted.values} if the \code{tp}'s are part of the
    fitted subset \code{result$control$subset}.
  }
  \item{which.start}{
    Which initial parameter values should be used when successively
    refitting the model to subsets of the data (up to time point
    \code{tp[1]}, up to \code{tp[1]+1}, ...) if \code{type="rolling"}?
    Default (\code{"current"}) is to use the parameter estimates from the
    previous time point, and \code{"final"} means to always
    use the estimates from \code{result} as initial values.
    Alternatively, \code{which.start} can be a list of \code{start}
    values as expected by \code{\link{hhh4}}, which then replace
    the corresponding estimates from \code{result} as initial values.
    This argument is ignored for \dQuote{non-rolling} \code{type}s.
  }
  \item{keep.estimates}{
    logical indicating if parameter estimates and log-likelihoods from
    the successive fits should be returned.
  }
  \item{verbose}{
    non-negative integer (usually in the range \code{0:3}) specifying
    the amount of tracing information to output.
    During \code{hhh4} model updates, the following verbosity is used:
    \code{0} if \code{cores > 1}, otherwise \code{verbose-1} if there
    is more than one time point to predict, otherwise \code{verbose}.
  }
  \item{cores}{the number of cores to use when computing
    the predictions for the set of time points \code{tp} in parallel
    (with \code{\link[parallel]{mclapply}}).
    Note that parallelization is not possible in the default setting
    \code{type="rolling"} and \code{which.start="current"} (use
    \code{which.start="final"} for this to work).}

  \item{object}{an object of class \code{"oneStepAhead"}.}
  \item{parm}{unused (argument of the generic).}
  \item{level}{required confidence level of the prediction interval.}
  \item{probs}{numeric vector of probabilities with values in [0,1].}
  \item{unit}{single integer or character selecting a unit for which to
    produce the plot.}
  \item{start}{
    x-coordinate of the first prediction. If \code{start=NULL}
    (default), this is derived from \code{x}.
  }
  \item{means.args}{
    if a list (of graphical parameters for \code{\link{lines}}),
    the point predictions (from \code{x$pred}) are added to the plot.
  }
  
  \item{x}{an object of class \code{"oneStepAhead"} or \code{"hhh4"}.}
  \item{which}{character vector determining which scores to compute.
    The package \pkg{surveillance} implements the following proper
    scoring rules: logarithmic score (\code{"logs"}), ranked probability
    score (\code{"rps"}), Dawid-Sebastiani score (\code{"dss"}), and
    squared error score (\code{"ses"}). The normalized SES
    (\code{"nses"}) is also available but it is improper and hence not
    computed by default.\cr
    It is possible to name own scoring rules in \code{which}. These
    must be functions of \code{(x, mu, size)}, vectorized in all arguments
    (time x unit matrices) except that \code{size} is \code{NULL}
    in case of a Poisson model.
    See the available scoring rules for guidance, e.g., \code{\link{dss}}.
  }
  \item{subset}{
    subset of time points for which to calculate the scores
    (or test calibration, or produce the PIT histogram, respectively).
    Defaults to the subset used for fitting the model.}
  \item{units}{integer or character vector indexing the units for which
    to compute the scores (or the calibration test or the PIT histogram,
    respectively). By default, all units are considered.}
  \item{sign}{logical indicating if the function should also return
    \code{sign(x-mu)}, i.e., the sign of the difference between
    the observed counts and corresponding predictions.
    This does not really make sense when averaging over multiple
    \code{units} with \code{individual=FALSE}.}
  \item{individual}{logical indicating if the individual scores of the
    \code{units} should be returned. By default (\code{FALSE}), the
    individual scores are averaged over all \code{units}.}
  \item{reverse}{logical indicating if the rows (time points) should be
    reversed in the result. The long-standing but awkward default was to
    do so for the \code{oneStepAhead}-method. This has changed in
    version 1.16.0, so time points are no longer reversed by default.}
  \item{\dots}{Unused by the \code{quantile}, \code{confint} and
    \code{scores} methods.\cr
    The \code{plot}-method passes further arguments to the
    \code{\link{fanplot}} function, e.g., \code{fan.args},
    \code{observed.args}, and \code{key.args} can be used to modify the
    plotting style.\cr
    For the \code{calibrationTest}-method, further arguments are passed
    to \code{\link{calibrationTest.default}}, e.g., \code{which} to
    select a scoring rule.\cr
    For the \code{pit}-methods, further arguments are passed to
    \code{\link{pit.default}}.}
}

\value{
\code{oneStepAhead} returns a list (of class \code{"oneStepAhead"})
with the following components:
\item{pred}{one-step-ahead predictions in a matrix, where each row
  corresponds to one of the time points requested via the argument
  \code{tp}, and which has \code{ncol(result$stsObj)}
  unit-specific columns. The rownames indicate the predicted time points
  and the column names are identical to \code{colnames(result$stsObj)}.}
\item{observed}{matrix with observed counts at the predicted time
  points. It has the same dimensions and names as \code{pred}.}
\item{psi}{in case of a negative-binomial model, a matrix of the
  estimated overdispersion parameter(s) at each time point on 
  the internal -log-scale (1 column if \code{"NegBin1"},
  \code{ncol(observed)} columns if \code{"NegBinM"} or shared overdispersion). 
  For a \code{"Poisson"} model, this component is \code{NULL}.}
\item{allConverged}{logical indicating if all successive fits
  converged.}

If \code{keep.estimates=TRUE}, there are the following additional elements:
\item{coefficients}{matrix of estimated regression parameters from the successive fits.}
\item{Sigma.orig}{matrix of estimated variance parameters from the successive fits.}
\item{logliks}{matrix with columns \code{"loglikelihood"} and
  \code{"margll"} with their obvious meanings.}

The \code{quantile}-method computes quantiles of the one-step-ahead
forecasts. If there is only one unit, it returns a tp x prob matrix,
otherwise a tp x unit x prob array.
The \code{confint}-method is a convenient wrapper with \code{probs} set
according to the required confidence level.

The function \code{scores} computes the scoring rules specified in the
argument \code{which}.
If multiple \code{units} are selected and \code{individual=TRUE}, the
result is an array of dimensions
\code{c(nrow(pred),length(units),5+sign)} (up to \pkg{surveillance}
1.8-0, the first two dimensions were collapsed to give a matrix).
Otherwise, the result is a matrix with \code{nrow(pred)} rows and
\code{5+sign} columns. If there is only one predicted time point, the
first dimension is dropped in both cases.


The \code{\link{calibrationTest}}- and \code{\link{pit}}-methods are
just convenient wrappers around the respective default methods.
}

\references{
  Czado, C., Gneiting, T. and Held, L. (2009):
    Predictive model assessment for count data.
    \emph{Biometrics}, \bold{65} (4), 1254-1261.
    \doi{10.1111/j.1541-0420.2009.01191.x}

  Paul, M. and Held, L. (2011):
    Predictive assessment of a non-linear random effects model for
    multivariate time series of infectious disease counts.
    \emph{Statistics in Medicine}, \bold{30} (10), 1118-1136.
    \doi{10.1002/sim.4177}
}

\author{
Sebastian Meyer and Michaela Paul
}

\seealso{
\code{vignette("hhh4")} and \code{vignette("hhh4_spacetime")}
}

\examples{
### univariate salmonella agona count time series

data("salmonella.agona")
## convert from old "disProg" to new "sts" class
salmonella <- disProg2sts(salmonella.agona)

## generate formula for temporal and seasonal trends
f.end <- addSeason2formula(~1 + t, S=1, period=52)
model <- list(ar = list(f = ~1), end = list(f = f.end), family = "NegBin1")
## fit the model
result <- hhh4(salmonella, model)

## do sequential one-step-ahead predictions for the last 5 weeks
pred <- oneStepAhead(result, nrow(salmonella)-5, type="rolling",
                     which.start="final", verbose=FALSE)
pred
quantile(pred)
confint(pred)

## simple plot of the 95% one-week-ahead prediction interval
## and point forecasts
plot(pred, probs = c(2.5,97.5)/100, means.args = list())

\dontshow{
## test equivalence of parallelized version
if (.Platform$OS.type == "unix" && isTRUE(parallel::detectCores() > 1))
    stopifnot(identical(pred,
        oneStepAhead(result, nrow(salmonella)-5, type="rolling",
                     which.start="final", verbose=FALSE, cores=2)))
}

## note: oneStepAhead(..., type="final") just means fitted values
stopifnot(identical(
    unname(oneStepAhead(result, nrow(salmonella)-5,
                        type="final", verbose=FALSE)$pred),
    unname(tail(fitted(result), 5))))


## compute scores of the one-step-ahead predictions
(sc <- scores(pred))

## the above uses the scores-method for "oneStepAhead" predictions,
## which is a simple wrapper around the default method:
scores(x = pred$observed, mu = pred$pred, size = exp(pred$psi))

## scores with respect to the fitted values are similar
(scFitted <- scores(result, subset = nrow(salmonella)-(4:0)))

\dontshow{
## test that scFitted is equivalent to scores(oneStepAhead(..., type = "final"))
stopifnot(all.equal(
    scFitted,
    scores(oneStepAhead(result, nrow(salmonella)-5, type="final", verbose=FALSE)),
    check.attributes = FALSE))
}


## test if the one-step-ahead predictions are calibrated
calibrationTest(pred)  # p = 0.8746

## the above uses the calibrationTest-method for "oneStepAhead" predictions,
## which is a simple wrapper around the default method:
calibrationTest(x = pred$observed, mu = pred$pred, size = exp(pred$psi))

## we can also test calibration of the fitted values
## using the calibrationTest-method for "hhh4" fits
calibrationTest(result, subset = nrow(salmonella)-(4:0))


## plot a (non-randomized) PIT histogram for the predictions
pit(pred)

## the above uses the pit-method for "oneStepAhead" predictions,
## which is a simple wrapper around the default method:
pit(x = pred$observed, pdistr = "pnbinom", mu = pred$pred, size = exp(pred$psi))


### multivariate measles count time series
## (omitting oneStepAhead forecasts here to keep runtime low)

data("measlesWeserEms")

## fit a hhh4 model with random effects in the endemic component
measlesModel <- list(
    end = list(f = addSeason2formula(~0 + ri(type="iid"))),
    ar = list(f = ~1),
    ne = list(f = ~1, weights = W_powerlaw(maxlag = 5)),
    family = "NegBin1")
measlesFit <- hhh4(measlesWeserEms, control = measlesModel)

## plot fitted mean components
plot(measlesFit, units = NULL)

## assess overall (in-sample) calibration of the model, i.e.,
## if the observed counts are from the fitted NegBin distribution
calibrationTest(measlesFit) # default is DSS (not suitable for low counts)
calibrationTest(measlesFit, which = "rps")  # p = 0.8267
calibrationTest(measlesFit, which = "logs") # p = 0.636

## to assess calibration in the second year for a specific district
calibrationTest(measlesFit, subset = 53:104, units = "03452", which = "rps")
pit(measlesFit, subset = 53:104, units = "03452")


### For a more sophisticated multivariate analysis of
### areal time series of influenza counts - data("fluBYBW") -
### see the (computer-intensive) demo("fluBYBW") script:
demoscript <- system.file(file.path("demo", "fluBYBW.R"),
                          package = "surveillance")
demoscript
#file.show(demoscript)
}

\keyword{univar}
\keyword{htest}
\keyword{dplot}
\keyword{ts}
