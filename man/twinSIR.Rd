\encoding{latin1}
\name{twinSIR}
\alias{twinSIR}

\title{
Fit an Additive-Multiplicative Intensity Model for SIR Data
}

\description{
\code{twinSIR} is used to fit additive-multiplicative intensity models for
epidemics as described in \enc{Höhle}{Hoehle} (2009).  Estimation is driven 
by (penalized) maximum likelihood in the point process frame work.  Optimization 
(maximization) of the (penalized) likelihood function is performed by means of 
\code{\link{optim}}.
The implementation is illustrated in Meyer et al. (2017, Section 4),
see \code{vignette("twinSIR")}.
}

\usage{
twinSIR(formula, data, weights, subset,
        knots = NULL, nIntervals = 1, lambda.smooth = 0, penalty = 1,
        optim.args = list(), model = TRUE, keep.data = FALSE)
}

\arguments{
  \item{formula}{
    an object of class \code{"\link{formula}"} (or one that can be coerced to
    that class): a symbolic description of the intensity model to be estimated.  
    The details of the model specification are given below.
  }
  \item{data}{
    an object inheriting from class \code{"\link{epidata}"}.
  }
  \item{weights}{
    an optional vector of weights to be used in the fitting process.  Should be
    \code{NULL} (the default, i.e. all observations have unit weight) or a
    numeric vector.
  }
  \item{subset}{
    an optional vector specifying a subset of observations to be used in the
    fitting process.  The subset \code{atRiskY == 1} is automatically chosen,
    because the likelihood only depends on those observations.
  }
  \item{knots}{
    numeric vector or \code{NULL} (the default).  Specification of the knots,
    where we suppose a step of the log-baseline.  With the current 
    implementation, these must be existing \code{"stop"} time points in the
    selected \code{subset} of the \code{data}, which is always
    restricted to \code{atRiskY == 1} rows.
    The intervals of constant log-baseline 
    hazard rate then are \eqn{(minTime;knots_1]}, \eqn{(knots_1;knots_2]},
    \ldots, \eqn{(knots_K;maxTime]}.
    By default, the \code{knots} are automatically chosen at the quantiles of
    the infection time points such that \code{nIntervals} intervals result.
    Non-NULL \code{knots} take precedence over \code{nIntervals}.
  }
  \item{nIntervals}{
    the number of intervals of constant log-baseline hazard.  Defaults to 1,
    which means an overall constant log-baseline hazard will be fitted.
  }
  \item{lambda.smooth}{
    numeric, the smoothing parameter \eqn{\lambda}.  By default it is 0 which
    leads to unpenalized likelihood inference.
    In case \code{lambda.smooth=-1}, the automatic smoothing parameter
    selection based on a mixed model approach is used (cf.
    \enc{Höhle}{Hoehle}, 2009).
  }
  \item{penalty}{
    either a single number denoting the order of the difference used to penalize
    the log-baseline coefficients (defaults to 1), or a more specific penalty
    matrix \eqn{K} for the parameter sub-vector \eqn{\beta}. In case of
    non-equidistant knots -- usually the case when using quantile based
    knot locations -- only a 1st order differences penalty matrix as in
    Fahrmeir and Lang (2001) is implemented.
  }
  \item{optim.args}{
    a list with arguments passed to the \code{\link{optim}} function.
    Especially useful are the following ones:
    \describe{
    \item{\code{par}:}{
      to specify initial parameter values.  Those must be in the order
      \code{c(alpha, h0, beta)}, i.e. first the coefficients of the epidemic
      covariates in the same order as they appear in the \code{formula}, then
      the log-baseline levels in chronological order and finally the
      coefficients of the endemic covariates in the same order
      as they appear in the \code{cox} terms of the \code{formula}.  The default
      is to start with 1's for \code{alpha} and 0's for \code{h0} and
      \code{beta}.
    }
    \item{\code{control}:}{
      for more detailed \code{trace}-ing (default: 1), another \code{REPORT}-ing
      frequency if \code{trace} is positive (default: 10), higher \code{maxit}
      (maximum number of iterations, default: 300) or another \code{factr} value
      (default: 1e7, a lower value means higher precision).
    }
    \item{\code{method}:}{
      the optimization algorithm defaults to \code{"L-BFGS-B"} (for
      box-constrained optimization), if there are any epidemic (non-\code{cox})
      variables in the model, and to \code{"BFGS"} otherwise.
    }
    \item{\code{lower}:}{
      if \code{method = "L-BFGS-B"} this defines the lower bounds for the
      model coefficients.  By default, all effects \eqn{\alpha} of epidemic
      variables are restricted to be non-negative.  Normally, this is exactly
      what one would like to have, but there might be reasons for other lower
      bounds, see the Note below.
    }
    \item{\code{hessian}:}{
      An estimation of the Expected Fisher Information matrix is always 
      part of the return value of the function.  It might be interesting to see 
      the Observed Fisher Information (= negative Hessian at the maximum), too. 
      This will be additionally returned if \code{hessian = TRUE}.
    }
    }
  }
  \item{model}{
    logical indicating if the model frame, the \code{weights},
    \code{lambda.smooth}, the penalty matrix \eqn{K} and the list of used
    distance functions \code{f} (from \code{attributes(data)}) should be
    returned for further computation.  This defaults to \code{TRUE} as this
    information is necessary e.g. in the \code{profile} and \code{plot} 
    methods.
  }
  \item{keep.data}{
    logical indicating if the \code{"epidata"} object (\code{data})
    should be part of the return value. This is only necessary for use of the
    \code{\link[=simulate.twinSIR]{simulate}}-method for \code{"twinSIR"}
    objects.  The reason is that the \code{twinSIR} function only uses and
    stores the rows with \code{atRiskY == 1} in the \code{model} component, but
    for the simulation of new epidemic data one needs the whole data set with
    all individuals in every time block.  The default value is \code{FALSE}, so
    if you intent to use \code{simulate.twinSIR}, you have to set this to
    \code{TRUE}.
  }
}

\details{  
  A model is specified through the \code{formula}, which has the form
  
  \code{~ epidemicTerm1 + epidemicTerm2 + cox(endemicVar1) *
    cox(endemicVar2)},
  
  i.e. the right hand side has the usual form as in \code{\link{lm}} with
  some variables marked as being endemic by the special function
  \code{\link{cox}}.  The left hand side of the formula is empty and will be
  set internally to \code{cbind(start, stop, event)}, which is similar to
  \code{Surv(start, stop, event, type="counting")} in package \pkg{survival}.
  
  Basically, the additive-multiplicative model for the infection intensity
  \eqn{\lambda_i(t)} for individual \eqn{i} is
  \deqn{\lambda_i(t) = Y_i(t) * (e_i(t) + h_i(t))}
  where
  \describe{
    \item{Y\_i(t)}{
      is the at-risk indicator, indicating if individual \eqn{i} is
      \dQuote{at risk} of becoming infected at time point \eqn{t}.
      This variable is part of the event history \code{data}.
    }
    \item{e\_i(t)}{
      is the epidemic component of the infection intensity, defined as
      \deqn{e_i(t) = \sum_{j \in I(t)} f(||s_i - s_j||)}
      where \eqn{I(t)} is the set of infectious individuals just before time
      point \eqn{t}, \eqn{s_i} is the coordinate vector of individual \eqn{i}
      and the function \eqn{f} is defined as
      \deqn{f(u) = \sum_{m=1}^p \alpha_m B_m(u)}
      with unknown transmission parameters \eqn{\alpha} and known distance
      functions \eqn{B_m}. This set of distance functions results in the set of
      epidemic variables normally calculated by the converter function
      \code{\link{as.epidata}}, considering the equality
      \deqn{e_i(t) = \sum_{m=1}^p \alpha_m x_{im}(t)}
      with \eqn{x_{im}(t) = \sum_{j \in I(t)} B_m(||s_i - s_j||)} being the
      \eqn{m}'th epidemic variable for individual \eqn{i}.
    }
    \item{h\_i(t)}{
      is the endemic (\code{cox}) component of the infection intensity, defined
      as
      \deqn{h_i(t) = \exp(h_0(t) + z_i(t)' \beta)}
      where \eqn{h_0(t)} is the log-baseline hazard function, \eqn{z_i(t)}
      is the vector of endemic covariates of individual \eqn{i} and \eqn{\beta}
      is the vector of unknown coefficients.
      To fit the model, the log-baseline hazard function is approximated by a
      piecewise constant function with known knots, but unknown levels,
      which will be estimated. The approximation is specified by the arguments
      \code{knots} or \code{nIntervals}.
    }
  }
  
  If a big number of \code{knots} (or \code{nIntervals}) is chosen, the
  corresponding log-baseline parameters can be rendered identifiable by
  the use of penalized likelihood inference.  At present, it is the job
  of the user to choose an adequate value of the smoothing parameter
  \code{lambda.smooth}. Alternatively, a data driven
  \code{lambda.smooth} smoothing parameter selection based on a mixed
  model representation of an equivalent truncated power spline is offered (see
  reference for further details). The following two steps are iterated
  until convergence:
  \enumerate{
    \item Given fixed smoothing parameter, the penalized
    likelihood is optimized for the regression components using a L-BFGS-B
    approach
    \item Given fixed regression parameters, a Laplace approximation of the
    marginal likelihood for the smoothing parameter is numerically
    optimized.  
  }
  Depending on the data, convergence might take a couple of iterations.

  
  Note also that it is unwise to include endemic covariates with huge values,
  as they affect the intensities on the exponential scale (after
  multiplication by the parameter vector \eqn{\beta}).
  With large covariate values, the
  \code{optim} method "L-BFGS-B" will likely terminate due to an infinite 
  log-likelihood or score function in some iteration.
}

\value{
  \code{twinSIR} returns an object of class
  \code{"twinSIR"}, which is a list containing the following components:
  
  \item{coefficients}{a named vector of coefficients.}
  \item{loglik}{the maximum of the (penalized) log-likelihood function.}
  \item{counts}{the number of log-likelihood and score function evaluations.}
  \item{converged}{logical indicating convergence of the optimization
    algorithm.}
  \item{fisherinfo.observed}{if requested, the negative Hessian from
    \code{optim}.}
  \item{fisherinfo}{an estimation of the Expected Fisher Information matrix.}
  \item{method}{the optimization algorithm used.}
  \item{intervals}{a numeric vector (\code{c(minTime, knots, maxTime)})
    representing the consecutive intervals of constant log-baseline.}
  \item{nEvents}{a numeric vector containing the number of infections in each of
    the above \code{intervals}.}
  \item{model}{if requested, the model information used. This is a list with
    components \code{"survs"} (data.frame with the id, start, stop and event
    columns), \code{"X"} (matrix of the epidemic variables), \code{"Z"} (matrix
    of the endemic variables), \code{"weights"} (the specified \code{weights}), 
    \code{"lambda.smooth"} (the specified \code{lambda.smooth}), \code{"K"}
    (the penalty matrix used), and \code{"f"} and \code{"w"}
    (the functions to generate the used epidemic covariates).
    Be aware that the model only contains those rows with \code{atRiskY == 1}!}
  \item{data}{if requested, the supplied \code{"epidata"} \code{data}.}
  \item{call}{the matched call.}
  \item{formula}{the specified \code{formula}.}
  \item{terms}{the \code{terms} object used.}
}

\references{
  \enc{Höhle}{Hoehle}, M. (2009),
  Additive-multiplicative regression models for spatio-temporal
  epidemics, \emph{Biometrical Journal}, \bold{51} (6), 961-978.

  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
}

\author{
Michael \enc{Höhle}{Hoehle} and Sebastian Meyer
}

\note{
There are some restrictions to modelling the infection intensity
without a baseline hazard rate, i.e. without an intercept in the
\code{formula}.  
Reason: At some point, the optimization algorithm L-BFGS-B tries to set all 
transmission parameters \eqn{\alpha} to the boundary value 0 and to calculate
the (penalized) score function with this set of parameters (all 0).  The problem
then is that the values of the infection intensities \eqn{lambda_i(t)} are 0
for all \eqn{i} and \eqn{t} and especially at observed event times, which is 
impossible.  Without a baseline, it is not allowed to have all alpha's set to 0, 
because then we would not observe any infections.  Unfortunately, L-BFGS-B can 
not consider this restriction.  Thus, if one wants to fit a model without 
baseline hazard, the control parameter \code{lower} must be specified in 
\code{optim.args} so that some alpha is strictly positive, e.g.
\code{optim.args = list(lower = c(0,0.001,0.001,0))} and the initial parameter
vector \code{par} must not be the zero vector.
}

\seealso{
\code{\link{as.epidata}} for the necessary data input structure,
\code{\link{plot.twinSIR}} for plotting the path of the infection intensity,
\code{\link{profile.twinSIR}} for profile likelihood estimation.
and \code{\link{simulate.twinSIR}} for the simulation of epidemics following
the fitted model.

Furthermore, the standard extraction methods
\code{\link[=vcov.twinSIR]{vcov}}, \code{\link[=logLik.twinSIR]{logLik}},
\code{\link[=AIC.twinSIR]{AIC}} and
\code{\link[=extractAIC.twinSIR]{extractAIC}} are implemented for
objects of class \code{"twinSIR"}.
}

\examples{
data("hagelloch")
summary(hagelloch)

# simple model with an overall constant baseline hazard rate
fit1 <- twinSIR(~ household + cox(AGE), data = hagelloch)
fit1
summary(fit1)   # see also help("summary.twinSIR")
plot(fit1)      # see also help("plot.twinSIR")
checkResidualProcess(fit1)   # could be better

# fit a piecewise constant baseline hazard rate with 3 intervals using 
# _un_penalized ML and estimated coefs from fit1 as starting values 
fit2 <- twinSIR(~ household, data = hagelloch, nIntervals = 3,
                optim.args = list(par = coef(fit1)[c(1,2,2,2)]))
summary(fit2)

# fit a piecewise constant baseline hazard rate with 7 intervals
# using _penalized_ ML
fit3 <- twinSIR(~ household, data = hagelloch, nIntervals = 7,
                lambda.smooth = 0.1, penalty = 1)
summary(fit3)
checkResidualProcess(fit3)

# plot the estimated log-baseline levels
plot(x=fit2$intervals, y=coef(fit2)[c(2,2:4)], type="S", ylim=c(-6, -1))
lines(x=fit3$intervals, y=coef(fit3)[c(2,2:8)], type="S", col=2)
legend("right", legend=c("unpenalized 3", "penalized 7"), lty=1, col=1:2, bty="n")


## special use case: fit the model to a subset of the events only,
## while preserving epidemic contributions from the remainder
## (maybe some buffer area nodes)
fit_subset <- twinSIR(~ household, data = hagelloch, subset = CL=="preschool")
summary(fit_subset)

\dontshow{
    ## the eventTimes attribute was wrong in surveillance <= 1.15.0
    stopifnot(
        length(residuals(fit_subset)) == sum(fit_subset$model$survs$event)
    )
}
}

\keyword{models}
\keyword{optimize}
