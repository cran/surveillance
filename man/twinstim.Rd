\encoding{latin1}
\name{twinstim}
\alias{twinstim}

\title{
  Fit a Two-Component Spatio-Temporal Point Process Model
}

\description{
  A \code{twinstim} model as described in Meyer et al. (2012) is fitted to
  marked spatio-temporal point process data. This constitutes a
  regression approach for conditional intensity function modelling.
  The implementation is illustrated in Meyer et al. (2017, Section 3),
  see \code{vignette("twinstim")}.
}

\usage{
twinstim(endemic, epidemic, siaf, tiaf, qmatrix = data$qmatrix, data,
         subset, t0 = data$stgrid$start[1], T = tail(data$stgrid$stop,1),
         na.action = na.fail, start = NULL, partial = FALSE,
         epilink = "log", control.siaf = list(F = list(), Deriv = list()),
         optim.args = list(), finetune = FALSE,
         model = FALSE, cumCIF = FALSE, cumCIF.pb = interactive(),
         cores = 1, verbose = TRUE)
}

\arguments{
  \item{endemic}{
    right-hand side formula for the exponential (Cox-like
    multiplicative) endemic component. May contain offsets (to be marked
    by the special function \code{offset}).  If omitted or \code{~0}
    there will be no endemic component in the model.  A type-specific
    endemic intercept can be requested by including the term
    \code{(1|type)} in the formula.
  }
  \item{epidemic}{
    formula representing the epidemic model for the event-specific
    covariates (marks) determining infectivity. Offsets are not
    implemented here. If omitted or \code{~0} there will be no epidemic
    component in the model.
  }
  \item{siaf}{
    spatial interaction function. Possible specifications are:
    \itemize{
      \item \code{NULL} or missing, corresponding to
      \code{siaf.constant()}, i.e. spatially homogeneous 
      infectivity independent of the distance from the host
      \item a list as returned by \code{\link{siaf}} or, more commonly,
      generated by a predefined interaction function such as
      \code{\link{siaf.gaussian}} as in Meyer et al. (2012) or
      \code{\link{siaf.powerlaw}} as in Meyer and Held (2014).
      The latter requires unique event locations, possibly after random
      tie-breaking (\code{\link{untie}}) or imputation of
      interval-censored locations.
      \code{\link{siaf.exponential}} is a simpler alternative.
      \item a numeric vector corresponding to the knots of a step
      function, i.e. the same as \code{\link{siaf.step}(knots)}
    }
    If you run into \dQuote{false convergence} with a non-constant
    \code{siaf} specification, the numerical accuracy of the cubature
    methods is most likely too low (see the \code{control.siaf}
    argument).
  }
  \item{tiaf}{
    temporal interaction function. Possible specifications are:
    \itemize{
      \item \code{NULL} or missing, corresponding to
      \code{tiaf.constant()}, i.e. time-constant infectivity
      \item a list as returned by \code{\link{tiaf}} or by a
      predefined interaction function such as
      \code{\link{tiaf.exponential}}
      \item a numeric vector corresponding to the knots of a step
      function, i.e. the same as \code{\link{tiaf.step}(knots)}
    }
  }
  \item{qmatrix}{
    square indicator matrix (0/1 or \code{FALSE}/\code{TRUE}) for
    possible transmission between the event types. The matrix will be
    internally converted to \code{logical}. Defaults to the \eqn{Q} matrix
    specified in \code{data}.
  }
  \item{data}{
    an object of class \code{"\link{epidataCS}"}.
  }
  \item{subset}{
    an optional vector evaluating to logical indicating a subset of
    \code{data$events} to keep. Missing values are
    taken as \code{FALSE}. The expression is evaluated in the context of the
    \code{data$events@data} \code{data.frame}, i.e. columns of this
    \code{data.frame} may be referenced directly by name.
  }
  \item{t0, T}{
    events having occurred during (-Inf;t0] are regarded as part of the
    prehistory \eqn{H_0} of the process. Only events that occurred in the
    interval (t0; T] are considered in the likelihood.
    The time point \code{t0} (\code{T}) must
    be an element of \code{data$stgrid$start} (\code{data$stgrid$stop}).
    The default time range covers the whole spatio-temporal grid
    of endemic covariates.
  }
  \item{na.action}{
    how to deal with missing values in \code{data$events}? Do not use
    \code{\link{na.pass}}. Missing values in the spatio-temporal grid
    \code{data$stgrid} are not accepted.
  }
  \item{start}{
    a named vector of initial values for (a subset of) the parameters.
    The names must conform to the conventions of \code{twinstim} to be
    assigned to the correct model terms. For instance,
    \code{"h.(Intercept)"} = endemic intercept,
    \code{"h.I(start/365)"} = coefficient of a linear time
    trend in the endemic component, \code{"h.factorB"} =
    coefficient of the level B of the factor variable \code{factor} in
    the endemic predictor, \code{"e.(Intercept)"} = epidemic intercept,
    \code{"e.VAR"} = coefficient of the epidemic term \code{VAR},
    \code{"e.siaf.2"} = second \code{siaf} parameter,
    \code{"e.tiaf.1"} = first \code{tiaf} parameter.
    Elements which don't match any of the model parameters are ignored.

    Alternatively, \code{start} may also be a named list with elements
    \code{"endemic"} or \code{"h"}, \code{"epidemic"} or \code{"e"},
    \code{"siaf"} or \code{"e.siaf"}, and \code{"tiaf"} or \code{"e.tiaf"},
    each of which containing a named numeric vector with the term labels
    as names (i.e. without the prefix \code{"h."}, \code{"e."}, etc).
    Thus, \code{start=list(endemic=c("(Intercept)"=-10))} is equivalent
    to \code{start=c("h.(Intercept)"=-10)}.
  }
  \item{partial}{
    logical indicating if a partial likelihood similar to the approach
    by Diggle et al. (2010) should be used (default is \code{FALSE}).
    Note that the partial likelihood implementation is not well tested.
  }
  \item{epilink}{
    a character string determining the link function to be used for the
    \code{epidemic} linear predictor of event marks. By default, the
    log-link is used. The experimental alternative
    \code{epilink = "identity"} (for use by \code{\link{epitest}}) does
    not guarantee the force of infection to be positive. If this leads
    to a negative total intensity (endemic + epidemic), the point
    process is not well defined (the log-likelihood will be
    \code{\link{NaN}}).
  }
  \item{control.siaf}{
    a list with elements \code{"F"} and \code{"Deriv"}, which are lists
    of extra arguments passed to the functions \code{siaf$F} and
    \code{siaf$Deriv}, respectively.\cr
    These arguments control the accuracy of the cubature routines from
    package \pkg{polyCub} involved in non-constant \code{siaf}
    specifications, e.g., the bandwidth of the midpoint rule
    \code{\link[polyCub]{polyCub.midpoint}}, the number of Gaussian
    quadrature points for \code{\link[polyCub]{polyCub.SV}}, or the relative
    tolerance of \code{\link{integrate}} in \code{\link[polyCub]{polyCub.iso}}.\cr
    For instance, \code{\link{siaf.gaussian}(F.adaptive = TRUE)} uses the
    midpoint-cubature \code{\link[polyCub]{polyCub.midpoint}} with an
    adaptive bandwidth of \code{eps=adapt*sd} to numerically integrate the
    kernel \eqn{f(\bold{s})}, and the default \code{adapt} value (0.1)
    can be overwritten by setting \code{control.siaf$F$adapt}.
    However, the default version \code{siaf.gaussian()}
    as well as \code{\link{siaf.powerlaw}()} and friends use
    \code{\link[polyCub]{polyCub.iso}} and thus accept control arguments for the
    standard \code{\link{integrate}} routine (such as \code{rel.tol})
    via \code{control.siaf$F} and \code{control.siaf$Deriv}.\cr
    This argument list is ignored in the case
    \code{siaf=siaf.constant()} (which is the default if \code{siaf} is
    unspecified).
  }
  \item{optim.args}{
    an argument list passed to \code{\link{optim}}, or \code{NULL}, in
    which case no optimization will be performed
    but the necessary functions will be returned in a list (similar to
    what is returned if \code{model = TRUE}).

    Initial values for the parameters may be given as list element
    \code{par} in the order \code{(endemic, epidemic, siaf, tiaf)}.
    If no initial values are provided, crude estimates will be used for
    the endemic intercept and the Gaussian kernel, -9 for the epidemic
    intercept, and zeroes for the remaining parameters.
    Any initial values given in the \code{start} argument
    take precedence over those in \code{par}.
    
    Note that \code{optim} receives the negative log-likelihood for 
    minimization (thus, if used, \code{optim.args$control$fnscale} should be
    positive). The \code{hessian} argument defaults to \code{TRUE}, and
    in the \code{control} list, \code{trace}ing is enabled with
    \code{REPORT=1} by default. By setting
    \code{optim.args$control$trace = 0}, all output from the
    optimization routine is suppressed.

    For the \code{partial} likelihood, the analytic score function and
    the Fisher information are not implemented and the default is to use
    robust \code{method="Nelder-Mead"} optimization.
    
    There may be an extra component \code{fixed} in the
    \code{optim.args} list, which determines which parameters should stick
    to their initial values. This can be specified by a
    logical vector of the same length as the \code{par} component, by an
    integer vector indexing \code{par} or by a character vector following
    the \code{twinstim} naming conventions. Furthermore, if
    \code{isTRUE(fixed)}, then all parameters are fixed at their initial
    values and no optimization is performed.
    
    Importantly, the \code{method} argument in the \code{optim.args}
    list may also be \code{"nlminb"}, in 
    which case the \code{\link{nlminb}} optimizer is used. This is also the
    default for full likelihood inference.
    In this case, not only the score function but also the
    \emph{expected} Fisher information can be used during optimization (as
    estimated by what Martinussen and Scheike (2006, p. 64) call the
    \dQuote{optional variation process}, or see Rathbun (1996, equation
    (4.7))). In our experience this gives better convergence than
    \code{optim}'s methods.
    For \code{method="nlminb"}, the following parameters of the
    \code{optim.args$control} list may be named like for
    \code{optim} and are renamed appropriately:
    \code{maxit} (-> \code{iter.max}), \code{REPORT} (-> \code{trace},
    default: 1), \code{abstol} (-> \code{abs.tol}), and
    \code{reltol} (-> \code{rel.tol}, default: \code{1e-6}).
    For \code{nlminb}, a logical \code{hessian} argument (default:
    \code{TRUE}) indicates if the negative \emph{expected} Fisher
    information matrix should be used as the Hessian during optimization
    (otherwise a numerical approximation is used).

    Similarly, \code{method="nlm"} should also work but is not
    recommended here.    
  }
  \item{finetune}{
    logical indicating if a second maximisation should be performed with
    robust Nelder-Mead \code{optim} using the resulting
    parameters from the first maximisation as starting point. This
    argument is only considered if \code{partial = FALSE} and the
    default is to not conduct a second maximization (in most cases this
    does not improve upon the MLE).
  }
  \item{model}{
    logical indicating if the model environment should be kept with the
    result, which is required for
    \code{\link[=intensityplot.twinstim]{intensityplot}}s and
    \code{\link[=R0.twinstim]{R0}(..., trimmed = FALSE)}.
    Specifically, if \code{model=TRUE}, the return value will have the
    evaluation environment set as its \code{\link{environment}},
    and the returned \code{functions} element will contain the
    log-likelihood function (or partial log-likelihood function, if
    \code{partial = TRUE}), and optionally the score and the expected
    Fisher information functions (not for the partial likelihood, and
    only if \code{siaf} and \code{tiaf} provide the necessary
    derivatives).\cr
    Note that fitted objects with a model environment might consume
    quite a lot of memory since they contain the \code{data}.
  }
  \item{cumCIF}{
    logical (default: \code{FALSE}) indicating whether to
    calculate the fitted cumulative ground intensity at event times.
    This is the residual process, see \code{\link{residuals.twinstim}}.
  }
  \item{cumCIF.pb}{
    logical indicating if a progress bar should be shown
    during the calculation of \code{cumCIF}. Defaults to do so in an
    interactive \R session, and will be \code{FALSE} if \code{cores != 1}.
  }
  \item{cores}{
    number of processes to use in parallel operation. By default
    \code{twinstim} runs in single-CPU mode. Currently, only the
    \pkg{multicore}-type of parallel computing via forking is
    supported, which is not available on Windows, see
    \code{\link[parallel]{mclapply}} in package \pkg{parallel}.
    Note that for a \pkg{memoise}d \code{\link{siaf.step}} kernel,
    \code{cores=1} is fixed internally since parallelization would slow
    down model fitting significantly.
  }
  \item{verbose}{
    logical indicating if information should be printed during
    execution. Defaults to \code{TRUE}.
  }
}

\details{
  The function performs maximum likelihood inference
  for the additive-multiplicative spatio-temporal intensity model
  described in Meyer et al. (2012). It uses \code{\link{nlminb}} as the
  default optimizer and returns an object of class \code{"twinstim"}.
  Such objects have \code{print}, \code{\link[=plot.twinstim]{plot}} and
  \code{\link[=summary.twinstim]{summary}} methods.
  The \code{summary} output can be converted via corresponding
  \code{\link[=xtable.summary.twinstim]{xtable}} or
  \code{\link[=toLatex.summary.twinstim]{toLatex}} methods.
  Furthermore, the usual accessor methods are implemented, including
  \code{coef}, \code{vcov}, \code{logLik},
  \code{\link[=residuals.twinstim]{residuals}}, and
  \code{update}.
  Additional functionality is provided by the \code{\link{R0}} and
  \code{\link[=simulate.twinstim]{simulate}} methods.
}

\value{
  Returns an S3 object of class \code{"twinstim"}, which is a list with
  the following components:
  
  \item{coefficients}{vector containing the MLE.}
  \item{loglik}{value of the log-likelihood function at the MLE with a
    logical attribute \code{"partial"} indicating if the partial
    likelihood was used.}
  \item{counts}{number of log-likelihood and score evaluations during
    optimization.}
  \item{converged}{either \code{TRUE} (if the optimizer converged) or a
    character string containing a failure message.}
  \item{fisherinfo}{\emph{expected} Fisher information evaluated at the
    MLE. Only non-\code{NULL} for full likelihood inference
    (\code{partial = FALSE}) and if spatial and temporal interaction
    functions are provided with their derivatives.}
  \item{fisherinfo.observed}{observed Fisher information matrix
    evaluated at the value of the MLE. Obtained as the negative Hessian.
    Only non-\code{NULL} if \code{optim.args$method} is not
    \code{"nlminb"} and if it was requested by setting
    \code{hessian=TRUE} in \code{optim.args}.}
  \item{fitted}{fitted values of the conditional intensity function at the events.}
  \item{fittedComponents}{two-column matrix with columns \code{"h"} and
    \code{"e"} containing the fitted values of the endemic and epidemic
    components, respectively.\cr
    (Note that \code{rowSums(fittedComponents) == fitted}.)}
  \item{tau}{fitted cumulative ground intensities at the event times.
    Only non-\code{NULL} if \code{cumCIF = TRUE}.
    This is the \dQuote{residual process} of the model, see
    \code{\link{residuals.twinstim}}.}
  \item{R0}{estimated basic reproduction number for each event. This
    equals the spatio-temporal integral of the epidemic intensity over
    the observation domain (t0;T] x W for each event.}
  \item{npars}{vector describing the lengths of the 5 parameter
    subvectors: endemic intercept(s) \eqn{\beta_0(\kappa)}, endemic
    coefficients \eqn{\beta}, epidemic coefficients \eqn{\gamma},
    parameters of the \code{siaf} kernel, and parameters of the
    \code{tiaf} kernel.}
  \item{qmatrix}{the \code{qmatrix} associated with the epidemic
    \code{data} as supplied in the model call.}
  \item{bbox}{the bounding box of \code{data$W}.}
  \item{timeRange}{the time range used for fitting: \code{c(t0,T)}.}
  \item{formula}{a list containing the four main parts of the model
    specification: \code{endemic}, \code{epidemic}, \code{siaf}, and
    \code{tiaf}.}
  \item{xlevels}{a record of the levels of the factors used in fitting.}
  \item{control.siaf}{see the \dQuote{Arguments} section above.}
  \item{optim.args}{input optimizer arguments used to determine the MLE.}
  \item{functions}{if \code{model=TRUE} this is a \code{list} with
    components \code{ll}, \code{sc} and \code{fi}, which are functions
    evaluating the log-likelihood, the score function and the expected
    Fisher information for a parameter vector \eqn{\theta}. The
    \code{environment} of these function is the model environment, which
    is thus retained in the workspace if \code{model=TRUE}. Otherwise,
    the \code{functions} component is \code{NULL}.}
  \item{call}{the matched call.}
  \item{runtime}{the \code{\link{proc.time}}-queried time taken
    to fit the model, i.e., a named numeric vector of length 5 of class
    \code{"proc_time"}, with the number of \code{cores} set as
    additional attribute.}

  If \code{model=TRUE}, the model evaluation environment is assigned to
  this list and can thus be queried by calling \code{environment()} on
  the result.
}

\note{
  \code{twinstim} makes use of the \pkg{memoise} package if it is
  available -- and that is highly recommended for non-constant
  \code{siaf} specifications to speed up calculations. Specifically, the
  necessary numerical integrations of the spatial interaction function
  will be cached such that they are only calculated once for every
  state of the \code{siaf} parameters during optimization.
}

\references{
  Diggle, P. J., Kaimi, I. & Abellana, R. (2010):
  Partial-likelihood analysis of spatio-temporal point-process data.
  \emph{Biometrics}, \bold{66}, 347-354.

  Martinussen, T. and Scheike, T. H. (2006):
  Dynamic Regression Models for Survival Data.
  Springer.

  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{�}{ae}t
  M\enc{�}{ue}nchen.\cr
  Available as \url{https://epub.ub.uni-muenchen.de/11703/}
  
  Meyer, S., Elias, J. and H\enc{�}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}

  Meyer, S. and Held, L. (2014):
  Power-law models for infectious disease spread.
  \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.
  \doi{10.1214/14-AOAS743}

  Meyer, S., Held, L. and \enc{H�hle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}

  Rathbun, S. L. (1996):
  Asymptotic properties of the maximum likelihood estimator for
  spatio-temporal point processes.
  \emph{Journal of Statistical Planning and Inference}, \bold{51}, 55-74.
}

\author{
  Sebastian Meyer

  Contributions to this documentation by
  Michael H\enc{�}{oe}hle and Mayeul Kauffmann.
}

\seealso{
  \code{vignette("twinstim")}.
  
  There is a \code{\link{simulate.twinstim}} method,
  which simulates the point process based on the fitted \code{twinstim}.

  A discrete-space alternative is offered by the \code{\link{twinSIR}}
  modelling framework.
}

\examples{
# Load invasive meningococcal disease data
data("imdepi")


### first, fit a simple endemic-only model

m_noepi <- twinstim(
    endemic = addSeason2formula(~ offset(log(popdensity)) + I(start/365-3.5),
                                S=1, period=365, timevar="start"),
    data = imdepi, subset = !is.na(agegrp)
)

## look at the model summary
summary(m_noepi)

## there is no evidence for a type-dependent endemic intercept (LR test)
m_noepi_type <- update(m_noepi, endemic = ~(1|type) + .)
pchisq(2*c(logLik(m_noepi_type)-logLik(m_noepi)), df=1, lower.tail=FALSE)


### add an epidemic component with just the intercept, i.e.
### assuming uniform dispersal in time and space up to a distance of
### eps.s = 200 km and eps.t = 30 days (see summary(imdepi))

m0 <- update(m_noepi, epidemic=~1, model=TRUE)

## summarize the model fit
summary(m0, correlation = TRUE, symbolic.cor = TRUE)

## the default confint-method can be used for Wald-CI's
confint(m0, level=0.95)

## same "untrimmed" R0 for every event (simple epidemic intercept model)
summary(R0(m0, trimmed=FALSE))

## plot the path of the fitted total intensity
plot(m0, "total intensity", tgrid=500)

if (surveillance.options("allExamples")) {
## extract "residual process" integrating over space (takes some seconds)
res <- residuals(m0)
# if the model describes the true CIF well _in the temporal dimension_,
# then this residual process should behave like a stationary Poisson
# process with intensity 1
plot(res, type="l"); abline(h=c(0, length(res)), lty=2)
# easier, with CI and serial correlation:
checkResidualProcess(m0)
}

\dontrun{
  ## NB: in contrast to nlminb(), optim's BFGS would miss the
  ##     likelihood maximum wrt the epidemic intercept
  m0_BFGS <- update(m_noepi, epidemic=~1, optim.args = list(method="BFGS"))
  format(cbind(nlminb=coef(m0), BFGS=coef(m0_BFGS)), digits=3, scientific=FALSE)
  m0_BFGS$fisherinfo   # singular Fisher information matrix here
  m0$fisherinfo
  logLik(m0_BFGS)
  logLik(m0)
  ## nlminb is more powerful since we make use of the analytical fisherinfo
  ## as estimated by the model during optimization, which optim cannot
}


### an epidemic-only model?
## for a purely epidemic model, all events must have potential source events
## (otherwise the intensity at the observed event would be 0)

## let's focus on the C-type for this example
imdepiC <- subset(imdepi, type == "C")
table(summary(imdepiC)$nSources)
## 106 events have no prior, close events (in terms of eps.s and eps.t)
try(twinstim(epidemic = ~1, data = imdepiC))  # detects this problem
## let's assume spatially unbounded interaction
imdepiC_infeps <- update(imdepiC, eps.s = Inf)
(s <- summary(imdepiC_infeps))
table(s$nSources)
## for 11 events, there is no prior event within eps.t = 30 days
## (which is certainly true for the first event)
plot(s$counter, main = "Number of infectious individuals over time (eps.t = 30)")
rug(imdepiC_infeps$events$time)
rug(imdepiC_infeps$events$time[s$nSources == 0], col = 2, lwd = 3)
## An endemic component would catch such events (from unobserved sources),
## otherwise a longer infectious period would need to be assumed and
## for the first event to happen, a prehistory is required (e.g., t0 = 31).
## As an example, we fit the data only until T = 638 (all events have ancestors)
m_epi <- twinstim(epidemic = ~1, data = imdepiC_infeps, t0 = 31, T = 638)
summary(m_epi)


if (surveillance.options("allExamples")) withAutoprint({

### full model with interaction functions (time-consuming)
## estimate an exponential temporal decay of infectivity
m1_tiaf <- update(m0, tiaf=tiaf.exponential())
plot(m1_tiaf, "tiaf", scaled=FALSE)

## estimate a step function for spatial interaction
summary(sourceDists <- getSourceDists(imdepi, "space"))
(knots <- quantile(sourceDists, c(5,10,20,40)/100))
m1_fstep <- update(m0, siaf=knots)
plot(m1_fstep, "siaf", scaled=FALSE)
rug(sourceDists, ticksize=0.02)

## estimate a continuously decreasing spatial interaction function,
## here we use the kernel of an isotropic bivariate Gaussian
m1 <- update(m0, siaf = siaf.gaussian())
AIC(m_noepi, m0, m1_fstep, m1)
summary(m1)  # e.siaf.1 is log(sigma), no test for H0: log(sigma) = 0
exp(confint(m1, "e.siaf.1"))  # a confidence interval for sigma
plot(m1, "siaf", scaled=FALSE)
## alternative: siaf.powerlaw() with eps.s=Inf and untie()d data,
##              see vignette("twinstim")

## add epidemic covariates
m2 <- update(m1, epidemic = ~ 1 + type + agegrp)
AIC(m1, m2)   # further improvement
summary(m2)
  
## look at estimated R0 values by event type
tapply(R0(m2), imdepi$events@data[names(R0(m2)), "type"], summary)

})
}

\keyword{models}
\keyword{optimize}
