\name{hhh4_simulate}
\alias{simulate.hhh4}

\title{Simulate \code{"hhh4"} Count Time Series}

\description{
   Simulates a multivariate time series of counts based on the
   Poisson/Negative Binomial model as described in Paul and Held (2011).
}

\usage{
\method{simulate}{hhh4}(object, nsim = 1, seed = NULL, y.start = NULL,
         subset = 1:nrow(object$stsObj), coefs = coef(object),
         components = c("ar","ne","end"), simplify = nsim>1, ...)
}

\arguments{
  \item{object}{
    an object of class \code{"\link{hhh4}"}.
  }
  \item{nsim}{
    number of time series to simulate. Defaults to \code{1}.
  }
  \item{seed}{
    an object specifying how the random number generator should be
    initialized for simulation (via \code{\link{set.seed}}). The
    initial state will also be stored as an attribute \code{"seed"} of
    the result. The original state of the \code{\link{.Random.seed}}
    will be restored at the end of the simulation.
    By default (\code{NULL}), neither initialization nor recovery will
    be done.
    This behaviour is copied from the \code{\link{simulate}.lm} method.
  }
  \item{y.start}{
    vector or matrix (with \code{ncol(object$stsObj)} columns) with
    starting counts for the epidemic components. 
    If \code{NULL}, the observed means in the respective units of the
    data in \code{object} during \code{subset} are used.
  }
  \item{subset}{
    time period in which to simulate data. Defaults to (and cannot
    exceed) the whole period defined by the underlying \code{"sts"}
    object.
  }
  \item{coefs}{
    coefficients used for simulation from the model in \code{object}.
    Default is to use the fitted parameters.
    Note that the \code{coefs}-vector must be in the same order and
    scaling as \code{coef(object)}, which especially means
    \code{reparamPsi = TRUE} (as per default when using the
    \code{coef}-method to extract the parameters).
    The overdispersion parameter in \code{coefs} is the inverse of the
    dispersion parameter \code{size} in \code{\link{rnbinom}}.
  }
  \item{components}{
    character vector indicating which components of the fitted model
    \code{object} should be active during simulation. For instance,
    a simulation with \code{components="end"} is solely based on the
    fitted endemic mean.
  }
  \item{simplify}{
    logical indicating if only the simulated counts (\code{TRUE}) or the
    full \code{"\linkS4class{sts}"} object (\code{FALSE}) should be
    returned for every replicate.
    By default a full \code{"sts"} object is returned iff \code{nsim=1}.
  }
  \item{\dots}{unused (argument of the generic).}
}

\details{
  Simulates data from a Poisson or a Negative Binomial model
  with mean
  \deqn{\mu_{it} = \lambda_{it} y_{i,t-1} + 
                   \phi_{it} \sum_{j \neq i} w_{ji} y_{j,t-1} + 
                   \nu_{it}}{%
        \mu_it = \lambda_it y_i,t-1 + 
                  \phi_it \sum_j w_ji y_j,t-1 + 
                  \nu_it}
  where
  \eqn{\lambda_{it}>0}, \eqn{\phi_{it}>0}, and \eqn{\nu_{it}>0} are 
  parameters which are modelled parametrically.
  The function uses the model and parameter estimates of the fitted
  \code{object} to simulate the time series.

  With the argument \code{coefs} it is possible to simulate from  
  the model as specified in \code{object}, but with different 
  parameter values.
}

\value{
  If \code{simplify=FALSE}: an object of class
  \code{"\linkS4class{sts}"} (\code{nsim = 1}) or a list of those
  (\code{nsim > 1}).

  If \code{simplify=TRUE}: an object of class
  \code{"hhh4sims"}, which is an array of dimension
  \code{c(length(subset), ncol(object$stsObj), nsim)}.
  The originally observed counts during the simulation period,
  \code{object$stsObj[subset,]}, are attached for reference
  (used by the \code{plot}-methods) as an attribute \code{"stsObserved"},
  and the initial condition \code{y.start} as attribute \code{"initial"}.
  The \code{[}-method for \code{"hhh4sims"} takes care of subsetting
  these attributes appropriately.
}

\references{
  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136
}
\author{
Michaela Paul and Sebastian Meyer
}

\seealso{
  \code{\link{plot.hhh4sims}} and \code{\link{scores.hhh4sims}}
  and the examples therein for \code{nsim > 1}.
}

\examples{
data(influMen)
# convert to sts class and extract meningococcal disease time series
meningo <- disProg2sts(influMen)[,2]

# fit model
fit <- hhh4(meningo, control = list(
              ar = list(f = ~ 1),
              end = list(f = addSeason2formula(~1, period = 52)),
              family = "NegBin1"))
plot(fit)

# simulate from model (generates an "sts" object)
simData <- simulate(fit, seed=1234)

# plot simulated data
plot(simData, main = "simulated data", xaxis.labelFormat=NULL)

# use simplify=TRUE to return an array of simulated counts
simCounts <- simulate(fit, seed=1234, simplify=TRUE)
dim(simCounts)  # nTime x nUnit x nsim
\dontshow{stopifnot(observed(simData) == c(simCounts))}
# plot the first year of simulated counts (+ initial + observed)
plot(simCounts[1:52,,], type = "time", xaxis.labelFormat = NULL)
# see help(plot.hhh4sims) for other plots, mainly useful for nsim > 1

# simulate from a Poisson instead of a NegBin model
# keeping all other parameters fixed at their original estimates
coefs <- replace(coef(fit), "overdisp", 0)
simData2 <- simulate(fit, seed=123, coefs = coefs)
plot(simData2, main = "simulated data: Poisson model", xaxis.labelFormat = NULL)

# simulate from a model with higher autoregressive parameter
coefs <- replace(coef(fit), "ar.1", log(0.9))
simData3 <- simulate(fit, seed=321, coefs = coefs)
plot(simData3, main = "simulated data: lambda = 0.5", xaxis.labelFormat = NULL)


## more sophisticated: simulate beyond initially observed time range

# extend data range by one year (non-observed domain), filling with NA values
nextend <- 52
timeslots <- c("observed", "state", "alarm", "upperbound", "populationFrac")
addrows <- function (mat, n) mat[c(seq_len(nrow(mat)), rep(NA, n)),,drop=FALSE]
extended <- Map(function (x) addrows(slot(meningo, x), n = nextend), x = timeslots)
# create new sts object with extended matrices
meningo2 <- do.call("sts", c(list(start = meningo@start, frequency = meningo@freq,
                                  map = meningo@map), extended))

# fit to the observed time range only, via the 'subset' argument
fit2 <- hhh4(meningo2, control = list(
              ar = list(f = ~ 1),
              end = list(f = addSeason2formula(~1, period = 52)),
              family = "NegBin1",
              subset = 2:(nrow(meningo2) - nextend)))
# the result is the same as before
stopifnot(all.equal(fit, fit2, ignore = c("stsObj", "control")))
\dontshow{
# one-week-ahead prediction only "works" for the first non-observed time point
# because the autoregressive component relies on non-missing past counts
oneStepAhead(fit2, tp = rep(nrow(meningo2)-nextend, 2), type = "final", verbose = FALSE)
# however, methods won't work as observed is NA
}
# long-term probabilistic forecast via simulation for non-observed time points
meningoSim <- simulate(fit2, nsim = 100, seed = 1,
                       subset = seq(nrow(meningo)+1, nrow(meningo2)),
                       y.start = tail(observed(meningo), 1))
apply(meningoSim, 1:2, function (ysim) quantile(ysim, c(0.1, 0.5, 0.9)))
# three plot types are available for "hhh4sims", see also ?plot.hhh4sims
plot(meningoSim, type = "time", average = median)
plot(meningoSim, type = "size", observed = FALSE)
if (requireNamespace("fanplot"))
    plot(meningoSim, type = "fan", means.args = list(),
         fan.args = list(ln = c(.1,.9), ln.col = 8))
}

\keyword{datagen}
