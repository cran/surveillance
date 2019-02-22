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
  \code{c(length(subset), ncol(object$stsObj), nsim)}, where the third
  dimension is dropped if \code{nsim=1} (yielding a matrix).
  The originally observed counts during the simulation period,
  \code{object$stsObj[subset,]}, are attached for reference
  (used by the \code{plot}-methods) as an attribute \code{"stsObserved"},
  and the initial condition \code{y.start} as attribute \code{"initial"}.
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
}

\examples{
data(influMen)
# convert to sts class and extract meningococcal disease time series
meningo <- disProg2sts(influMen)[,2]

# fit model
fit <- hhh4(meningo, control = list(ar = list(f = ~ 1),
            end = list(f = addSeason2formula(S = 1, period = 52)),
            family = "NegBin1"))
plot(fit)

# simulate from model
simData <- simulate(fit, seed=1234)

# plot simulated data
plot(simData, main = "simulated data", xaxis.labelFormat=NULL)

# consider a Poisson instead of a NegBin model
coefs <- coef(fit)
coefs["overdisp"] <- 0

simData2 <- simulate(fit, seed=123, coefs = coefs)
plot(simData2, main = "simulated data: Poisson model", xaxis.labelFormat = NULL)

# consider a model with higher autoregressive parameter
coefs <- coef(fit)
coefs[1] <- log(0.5)

simData3 <- simulate(fit, seed=321, coefs = coefs)
plot(simData3, main = "simulated data: lambda = 0.5", xaxis.labelFormat = NULL)
}

\keyword{datagen}
