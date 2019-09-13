\encoding{latin1}
\name{R0}
\alias{R0}
\alias{R0.twinstim}
\alias{R0.simEpidataCS}
\alias{simpleR0}

\title{Computes reproduction numbers from fitted models}

\description{
The S3 generic function \code{R0} defined in package \pkg{surveillance} is intended to
compute reproduction numbers from fitted epidemic models.
The package currently defines a method for the \code{"\link{twinstim}"} class, which
computes expected numbers of infections caused by infected individuals depending on the event type
and marks attached to the individual, which contribute to the infection pressure
in the epidemic predictor of that class.
There is also a method for simulated \code{"epidataCS"}
(just a wrapper for the \code{"twinstim"}-method).
}

\usage{
R0(object, ...)

\method{R0}{twinstim}(object, newevents, trimmed = TRUE, newcoef = NULL, ...)
\method{R0}{simEpidataCS}(object, trimmed = TRUE, ...)

simpleR0(object, eta = coef(object)[["e.(Intercept)"]],
         eps.s = NULL, eps.t = NULL, newcoef = NULL)
}

\arguments{
  \item{object}{A fitted epidemic model object for which an \code{R0} method exists.}
  \item{newevents}{
    an optional \code{data.frame} of events for which the reproduction
    numbers should be calculated.  If omitted, it is calculated for the
    original events from the fit.  In this case, if
    \code{trimmed = TRUE} (the default), the result is just
    \code{object$R0}; however, if \code{trimmed = FALSE}, the model
    environment is required, i.e. \code{object} must have been fitted
    with \code{model = TRUE}.
    
    For the \code{twinstim} method, \code{newevents} must at least
    contain the following columns:
    the event \code{time} (only for \code{trimmed = TRUE})
    and \code{type} (only for multi-type epidemics),
    the maximum interaction ranges \code{eps.t} and \code{eps.s},
    as well as columns for the marks and \code{stgrid} variables
    used in the epidemic component of the fitted \code{"twinstim"}
    \code{object} as stored in \code{formula(object)$epidemic}.
    For \code{trimmed} R0 values, \code{newevents} must additionally
    contain the components \code{.influenceRegion} and, if using the
    \code{Fcircle} trick in the \code{siaf} specification, also
    \code{.bdist} (cf. the hidden columns in the \code{events}
    component of class \code{"epidataCS"}).
  }
  \item{trimmed}{
    logical indicating if the individual reproduction numbers should be
    calculated by integrating the epidemic intensities over the
    observation period and region only (\code{trimmed = TRUE}) or over
    the whole time-space domain R+ x R^2 (\code{trimmed = FALSE}). By
    default, if \code{newevents} is missing, the trimmed \code{R0}
    values stored in \code{object} are returned. Trimming means that
    events near the (spatial or temporal) edges of the observation
    domain have lower reproduction numbers (ceteris paribus) because
    events outside the observation domain are not observed.
  }
  \item{newcoef}{
    the model parameters to use when calculating reproduction numbers.
    The default (\code{NULL}) is to use the MLE \code{coef(object)}.
    This argument mainly serves the construction of Monte Carlo
    confidence intervals by evaluating \code{R0} for parameter vectors
    sampled from the asymptotic multivariate normal distribution of the MLE,
    see Examples.
  }
  \item{\dots}{additional arguments passed to methods.
    Currently unused for the \code{twinstim} method.}

  \item{eta}{a value for the epidemic linear predictor, see details.}
  \item{eps.s,eps.t}{the spatial/temporal radius of interaction.
    If \code{NULL} (the default), the original value from the data is
    used if this is unique and an error is thrown otherwise.}
}

\details{
For the \code{"\link{twinstim}"} class, the individual-specific expected
number \eqn{\mu_j} of infections caused by individual (event) \eqn{j}
inside its theoretical (untrimmed) spatio-temporal range of interaction
given by its \code{eps.t} (\eqn{\epsilon}) and \code{eps.s}
(\eqn{\delta}) values is defined as follows (cf. Meyer et al, 2012):
\deqn{\mu_j = e^{\eta_j} \cdot
  \int_{b(\bold{0},\delta)} f(\bold{s}) d\bold{s} \cdot
  \int_0^\epsilon g(t) dt .}
Here, \eqn{b(\bold{0},\delta)} denotes the disc centred at (0,0)' with
radius \eqn{\delta}, \eqn{\eta_j} is the epidemic linear predictor,
\eqn{g(t)} is the temporal interaction function, and \eqn{f(\bold{s})}
is the spatial interaction function. For a type-specific
\code{twinstim}, there is an additional factor for the number of event
types which can be infected by the type of event \eqn{j} and the
interaction functions may be type-specific as well.

Alternatively to the equation above,
the \code{trimmed} (observed) reproduction numbers
are obtain by integrating over the observed infectious domains of the
individuals, i.e. integrate \eqn{f} over the intersection of the
influence region with the observation region \code{W}
(i.e. over \eqn{\{ W \cap b(\bold{s}_j,\delta) \} - \bold{s}_j})
and \eqn{g} over the intersection of the observed infectious period with
the observation period \eqn{(t_0;T]} (i.e. over
\eqn{(0; \min(T-t_j,\epsilon)]}).

The function \code{simpleR0} computes
\deqn{\exp(\eta) \cdot
  \int_{b(\bold{0},\delta)} f(\bold{s}) d\bold{s} \cdot 
  \int_0^{\epsilon} g(t) dt ,}
where \eqn{\eta} defaults to \eqn{\gamma_0} disregarding any epidemic
effects of types and marks. It is thus only
suitable for simple epidemic \code{\link{twinstim}} models with
\code{epidemic = ~1}, a diagonal (or secondary diagonal) \code{qmatrix},
and type-invariant interaction functions.
\code{simpleR0} mainly exists for use by \code{\link{epitest}}.

(Numerical) Integration is performed exactly as during the fitting of
\code{object}, for instance \code{object$control.siaf} is queried if
necessary.
}

\value{
  For the \code{R0} methods,
  a numeric vector of estimated reproduction numbers from the fitted
  model \code{object} corresponding to the rows of \code{newevents} (if
  supplied) or the original fitted events including events of the prehistory.

  For \code{simpleR0}, a single number (see details).
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}
}

\author{Sebastian Meyer}

\examples{
## load the 'imdepi' data and a model fit
data("imdepi", "imdepifit")

## calculate individual and type-specific reproduction numbers
R0s <- R0(imdepifit)
tapply(R0s, imdepi$events@data[names(R0s), "type"], summary)

## untrimmed R0 for specific event settings
refevent <- data.frame(agegrp = "[0,3)", type = "B", eps.s = Inf, eps.t = 30)
setting2 <- data.frame(agegrp = "[3,19)", type = "C", eps.s = Inf, eps.t = 14)
newevents <- rbind("ref" = refevent, "event2" = setting2)
(R0_examples <- R0(imdepifit, newevents = newevents, trimmed = FALSE))
stopifnot(all.equal(R0_examples[["ref"]],
                    simpleR0(imdepifit)))


### compute a Monte Carlo confidence interval

## use a simpler model with constant 'siaf' for speed
simplefit <- update(imdepifit, epidemic=~type, siaf=NULL, subset=NULL)

## we'd like to compute the mean R0's by event type
meanR0ByType <- function (newcoef) {
    R0events <- R0(simplefit, newcoef=newcoef)
    tapply(R0events, imdepi$events@data[names(R0events),"type"], mean)
}
(meansMLE <- meanR0ByType(newcoef=NULL))

## sample B times from asymptotic multivariate normal of the MLE
B <- 5  # CAVE: toy example! In practice this has to be much larger
set.seed(123)
parsamples <- MASS::mvrnorm(B, mu=coef(simplefit), Sigma=vcov(simplefit))

## for each sample compute the 'meanR0ByType'
meansMC <- apply(parsamples, 1, meanR0ByType)

## get the quantiles and print the result
cisMC <- apply(cbind(meansMLE, meansMC), 1, quantile, probs=c(0.025,0.975))
print(rbind(MLE=meansMLE, cisMC))


### R0 for a simple epidemic model
### without epidemic covariates, i.e., all individuals are equally infectious

mepi1 <- update(simplefit, epidemic = ~1, subset = type == "B",
                model = TRUE, verbose = FALSE)
## using the default spatial and temporal ranges of interaction
(R0B <- simpleR0(mepi1))  # eps.s=200, eps.t=30
stopifnot(identical(R0B, R0(mepi1, trimmed = FALSE)[[1]]))
## assuming smaller interaction ranges (but same infection intensity)
simpleR0(mepi1, eps.s = 50, eps.t = 15)
}

\keyword{methods}
\keyword{univar}
