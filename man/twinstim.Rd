\encoding{latin1}
\name{twinstim}
\alias{twinstim}

\title{
  Fit a Two-Component Spatio-Temporal Conditional Intensity Model
}

\description{
  A \code{twinstim} model as described in Meyer et al. (2012) is fitted to
  marked spatio-temporal point process data. This constitutes a
  regression approach for conditional intensity function modelling.
}

\usage{
twinstim(endemic, epidemic, siaf, tiaf, qmatrix = data$qmatrix, data,
         subset, t0 = data$stgrid$start[1], T = tail(data$stgrid$stop,1),
         na.action = na.fail, nCub, nCub.adaptive = FALSE, partial = FALSE,
         optim.args, finetune = FALSE,
         model = FALSE, cumCIF = TRUE, cumCIF.pb = TRUE)
}

\arguments{
  \item{endemic}{
    formula for the exponential (Cox-like multiplicative) endemic
    component. May contain offsets (to be marked by the special function
    \code{offset}).  If omitted or \code{~0} there will be no endemic
    component in the model.  A type-specific endemic intercept can be
    requested by including the term \code{(1|type)} in the formula.
  }
  \item{epidemic}{
    formula representing the epidemic model for the event-specific
    covariates (marks) determining infectivity. Offsets are not
    implemented here. If omitted or \code{~0} there will be no epidemic
    component in the model.
  }
  \item{siaf}{
    spatial interaction function. May be \code{NULL} or missing
    (corresponding to \code{siaf.constant()}, i.e. constant infectivity 
    independent of the spatial distance from the host), a list
    (specifying a continuous kernel, as e.g., returned by the functions
    \code{\link{siaf.gaussian}} or \code{\link{siaf.lomax}}), or numeric
    (knots of a step function).
  }
  \item{tiaf}{
    temporal interaction function. May be \code{NULL} (or missing,
    corresponding to constant infectivity independent of the temporal
    distance from the host), a list (specifying a continuous function,
    as e.g., returned by \code{\link{tiaf.exponential}}), or numeric
    (knots of a step function). 
  }
  \item{qmatrix}{
    square indicator matrix (0/1 or \code{TRUE}/\code{FALSE}) for
    possible transmission between the event types. Will be internally
    converted to logical. Defaults to the \eqn{Q} matrix specified in
    \code{data}.
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
    events having occured during (-Inf;t0] are regarded as part of the
    prehistory \eqn{H_0} of the process. Only events that occured in the
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
  \item{nCub, nCub.adaptive}{
    determine the accuracy of the numerical cubature for non-constant
    \code{siaf} specifications. These arguments are ignored in the case
    \code{siaf=siaf.constant()} (which is the default if \code{siaf} is
    unspecified).\cr
    If \code{siaf$Fcircle} and \code{siaf$effRange} are
    both specified and \code{nCub.adaptive = TRUE}, then \code{nCub} is
    interpreted as \code{siaf$effRange()}/\code{eps}, where \code{eps}
    is used as pixel width and height in the two-dimensional midpoint
    rule (see \code{\link{polyCub.midpoint}}). Thus \code{nCub} is the
    desired number of subdivions of the effective integration range in
    both dimensions (for the current \code{siaf} parameters during
    optimisation), and \code{eps} thus equals
    \code{siaf$effRange()}/\code{nCub}. This is intended and especially
    useful for the Gaussian \code{\link[=siaf.gaussian]{siaf}} kernel:
    if, e.g., 6\eqn{\sigma} is used as the effective integration range
    and \code{nCub = 24}, then a bandwidth of \eqn{\sigma / 4} is used
    in the midpoint rule.\cr 
    If \code{siaf$Fcircle} or \code{siaf$effRange} is unspecified or
    \code{nCub.adaptive = FALSE}, then \code{nCub} equals the above
    mentioned \code{eps} in \code{polyCub.midpoint} (in this case may
    also be a vector of length 2 to specify different bandwidths in
    the x and y directions, see \code{\link[spatstat]{as.mask}}.
  }
  \item{partial}{
    logical indicating if a partial log-likelihood similar to the approach
    by Diggle et al. (2010) should be used.
  }
  \item{optim.args}{
    \code{NULL} or an argument list passed to \code{\link{optim}} containing
    at least the element \code{par}, the start values of the parameters in
    the order \code{par = c(endemic, epidemic, siaf, tiaf)}.
    Note that \code{optim} receives the negative log-likelihood for minimization
    (thus, if used, \code{optim.args$control$fnscale} should be
    positive).

    There may be an extra component \code{fixed} in the
    \code{optim.args} list, which determines which parameters should stick
    to their initial (\code{par}) values. This can be specified by a
    logical vector of the same length as the \code{par} component, by an
    integer vector indexing \code{par} or by a character vector following
    the \code{twinstim} naming conventions (e.g. \code{"h.(Intercept)"} for
    the endemic intercept, \code{"e.siaf.1"} for the first \code{siaf}
    parameter).
    
    Importantly, the \code{method} argument in the \code{optim.args}
    list may also be \code{"nlminb"}, in 
    which case the \code{\link{nlminb}} optimizer is used. This is also the
    default. In this case, not only the score function but also the
    expected Fisher information is used during optimization (as
    estimated by what Martinussen and Scheike (2006, p. 64) call the
    \dQuote{optional variation process}, or see Rathbun (1996, equation
    (4.7))). In our experience this gives better convergence than \code{optim}.
    For \code{method="nlminb"}, the only parameters available in the
    \code{optim.args$control} list are
    \code{maxit}, \code{REPORT} (default: 5), \code{abstol}, and
    \code{reltol} (default: \code{1e-6}), which are passed appropriately
    to \code{\link{nlminb}}. There is no \code{hessian} argument for
    \code{nlminb}.
    
    If \code{optim.args} is \code{NULL}, then no optimization will be performed
    but the necessary functions will be returned in a list (similar to what is
    returned if \code{model = TRUE}).
  }
  \item{finetune}{
    logical indicating if a second maximisation should be performed with
    robust Nelder-Mead \code{optim} using the resulting
    parameters from the first maximisation as starting point.
  }
  \item{model}{
    logical. If \code{TRUE} the result contains an element \code{functions} with
    the log-likelihood function (or partial log-likelihood function, if
    \code{partial = TRUE}), and optionally the score function and
    the Fisher information function (not for the partial likelihood, and
    only if \code{siaf} and \code{tiaf} provide the necessary
    derivatives).
    The environment of those functions equals the evaluation environment
    of the fitting function, i.e. it is kept in the workspace and the
    necessary model frames are still available when \code{twinstim} has
    finished.  The environment is also set as an attribute of the
    \code{twinstim} result.
  }
  \item{cumCIF}{
    logical (default: \code{TRUE}) indicating whether to
    calculate the fitted cumulative ground intensity at event times.
    This is the residual process, see \code{\link{residuals.twinstim}}.
  }
  \item{cumCIF.pb}{
    logical indicating if a progress bar should be shown
    during the calculation of \code{cumCIF}. Defaults to \code{TRUE}.
  }
}

\details{
  The function performs maximum likelihood inference
  for the additive-multiplicative spatio-temporal intensity model
  described in Meyer et al. (2012). It uses 'nlminb' as the default optimizer
  (Newton-algorithm) and returns an object of class \code{twinstim}.
  Such ojects have working \code{print},
  \code{plot} and \code{summary} functions. The output of can again be
  processed by the \code{toLatex} function. Furthermore, the usual model
  fit functions such as \code{coef}, \code{vcov} and \code{logLik} work.
  A specific add on is the use of the function \code{R0}.
}

\value{
  Returns an S3 object of class \code{"twinstim"}, which is a list with
  the following components:
  
  \item{coefficients}{vector containing the MLE.}
  \item{loglik}{value of the log-likelihood function at the MLE with a
    logical attribute \code{"partial"} indicating if the partial
    likelihood was used.}
  \item{counts}{Number of log-likelihood and score evaluations during
    optimization.}
  \item{method}{Name of the optimizer used to determine the MLE.}
  \item{nCub, nCub.adaptive}{see the \dQuote{Arguments} section above.}
  \item{converged}{logical indicating if the optimizer converged.}
  \item{fisherinfo}{expected Fisher information evaluated at the
    MLE. Only part of the output for full likelihood inference
    (\code{partial = FALSE}) and if spatial and temporal interaction
    functions are provided with their derivatives.}
  \item{fisherinfo.observed}{observed Fisher information matrix
    evaluated at the value of the MLE. Obtained as the negative Hessian.
    This is only part of the result if \code{optim.args$method} is not
    \code{"nlminb"} (i.e. \code{optim} is used) and if it was requested by setting
    \code{hessian=TRUE} in \code{optim.args}.}
  \item{fitted}{fitted values of the conditional intensity function at the events.}
  \item{fittedComponents}{two-column matrix with columns \code{"h"} and
    \code{"e"} containing the fitted values of the endemic and epidemic
    components, respectively.\cr
    (Note that \code{rowSums(fittedComponents) == fitted}.)}
  \item{tau}{fitted cumulative ground intensities at the event times.
    Only calculated and returned if \code{cumCIF = TRUE}.
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
  \item{functions}{if \code{model=TRUE} this is a \code{list} with
    components \code{ll}, \code{sc} and \code{fi}, which are functions
    evaluating the log-likelihood, the score function and the expected
    Fisher information for a parameter vector \eqn{\theta}. The
    \code{environment} of these function is the model environment, which
    is thus retained in the workspace if \code{model=TRUE}.
  }
  \item{call}{the matched call.}
  \item{runtime}{contains the \code{\link{proc.time}} queried time taken
    to fit the model.}

  If \code{model=TRUE}, the model evaluation environment is assigned to
  this list and can thus be queried by calling \code{environment()} on
  the result.
}

\references{
  Diggle, P. J., Kaimi, I. & Abellana, R. (2010):
  Partial-likelihood analysis of spatio-temporal point-process data.
  \emph{Biometrics}, \bold{66}, 347-354.

  Martinussen, T. and Scheike, T. H. (2006):
  Dynamic Regression Models for Survival Data.
  Springer.
  
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}
  
  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{ä}{ae}t
  M\enc{ü}{ue}nchen.\cr
  Available as \url{http://epub.ub.uni-muenchen.de/11703/}

  Rathbun, S. L. (1996):
  Asymptotic properties of the maximum likelihood estimator for
  spatio-temporal point processes.
  \emph{Journal of Statistical Planning and Inference}, \bold{51}, 55-74.
}

\author{
  Sebastian Meyer with documentation contributions by Michael H\enc{ö}{oe}hle.
}

\seealso{
  \code{\link{twinSIR}} for a discrete space alternative.
  There is also a \code{\link{simulate.twinstim}} method,
  which simulates the point process based on the fitted \code{twinstim}.
}

\examples{
# Load invasive meningococcal disease data
data("imdepi")

### first, fit an endemic-only model

## As a start value for the endemic intercept use the crude estimate
## assuming the model only had a single-cell endemic component
## (rate of homogeneous Poisson process scaled for population density)
popdensity.overall <- with(subset(imdepi$stgrid, BLOCK == 1),
    weighted.mean(popdensity, area))
popdensity.overall   # pop. density in Germany is ~230 inhabitants/km^2
W.area <- with(subset(imdepi$stgrid, BLOCK == 1), sum(area))
W.area               # area of Germany is about 357100 km^2
# this should actually be the same as
sum(sapply(imdepi$W@polygons, slot, "area"))
# which here is not exactly the case because of polygon simplification

## start value for the endemic intercept
h.intercept <- with(summary(imdepi),
    log(nEvents/popdensity.overall/diff(timeRange)/W.area))

## fit the endemic-only model
m_noepi <- twinstim(
    endemic = ~ 1 + offset(log(popdensity)) + I(start/365) +
                sin(start * 2 * pi/365) + cos(start * 2 * pi/365),
    data = imdepi, subset = !is.na(agegrp),
    optim.args = list(par=c(h.intercept,rep(0,3)),
                      method="nlminb", control = list(REPORT=1)),
    model = FALSE, cumCIF = FALSE   # for reasons of speed
)

## look at the model summary
summary(m_noepi)

## type-dependent endemic intercept?
m_noepi_type <- update(m_noepi,
                       endemic = ~(1|type) + .,
                       optim.args=list(par=coef(m_noepi)[c(1,1:4)]))
summary(m_noepi_type)

# LR-test for a type-dependent intercept in the endemic-only model
pchisq(2*(logLik(m_noepi_type)-logLik(m_noepi)), df=1, lower.tail=FALSE)


### add an epidemic component with just the intercept, i.e.
### assuming uniform dispersal in time and space up to a distance of
### eps.s = 200 km and eps.t = 30 days (see summary(imdepi))

m0 <- update(m_noepi, epidemic = ~1,
             optim.args = list(par=c(coef(m_noepi),-10)), model = TRUE)

\dontrun{
## NOTE: optim's BFGS would not identify the epidemic intercept
m0_BFGS <- update(m_noepi, epidemic = ~1,
                  optim.args = list(par=c(coef(m_noepi),-10),
                                    method="BFGS"))
format(cbind(coef(m0), coef(m0_BFGS)), digits=3, scientific=FALSE)
m0_BFGS$fisherinfo   # singular Fisher information matrix here
m0$fisherinfo
logLik(m0_BFGS)
logLik(m0)
## nlminb is more powerful since it makes use of the fisherinfo as
## estimated by the model during optimization, which optim cannot
}

## summarize the model fit
s <- summary(m0, correlation = TRUE, symbolic.cor = TRUE)
s
# output the table of coefficients as LaTeX code
toLatex(s, digits=2, withAIC=FALSE)

## the default confint-method can be used for Wald-CI's
confint(m0, level=0.95)

## same (untrimmed) R0 estimate for every event (epidemic intercept only)
summary(R0(m0, trimmed=FALSE))

## plot the path of the total intensity or of the epidemic proportion
plot(m0, "total intensity", tgrid=1000)
plot(m0, "epidemic proportion", tgrid=1000)

\dontrun{
## extract "residual process" integrating over space (takes some seconds)
res <- residuals(m0)
# if the model describes the true CIF well _in the temporal dimension_,
# then this residual process should behave like a stationary Poisson
# process with intensity 1
plot(res, type="l"); abline(h=c(0, length(res)), lty=2)
# -> use the function checkResidualProcess
checkResidualProcess(m0)
}

## NB: the model fit environment is kept in the workspace with model=TRUE
sort(sapply(environment(m0), object.size))
# saving m0 to RData will include this model environment, thus might
# consume quiet a lot of disk space (as it does use memory)


### try with spatially decreasing interaction kernel
## Likelihood evaluation takes much longer than for constant spatial spread
## Here we use the kernel of an isotropic bivariate Gaussian with same
## standard deviation for both finetypes

m1 <- update(m0,
    siaf = siaf.gaussian(1, logsd = TRUE),
    optim.args = list(par=c(coef(m0), 2.8),  # exp(2.8) = sigma = 16 km
                      fixed="e.siaf.1"),
                      # for reasons of speed of the example, the siaf
                      # parameter log(sigma) is fixed to 2.8 here
    nCub = 24, nCub.adaptive = TRUE  # use adaptive eps=6*sigma/24=sigma/4
)

AIC(m_noepi, m0, m1)   # further improvement
summary(m1, test.iaf=FALSE)   # nonsense to test H0: log(sigma) = 0
plot(m1, "siaf", xlim=c(0,200), types=1)   # the same for types=2


### add epidemic covariates

m2 <- update(m1,
    epidemic = ~ 1 + type + agegrp,
    optim.args = list(par=c(coef(m0), 0, 0, 0, coef(m1)["e.siaf.1"]))
)

AIC(m1, m2)   # further improvement
summary(m2)

# look at estimated R0 values by event type
tapply(R0(m2), imdepi$events@data[names(R0(m2)), "type"], summary)
}

\keyword{models}
\keyword{optimize}
