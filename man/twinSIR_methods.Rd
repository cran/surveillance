\encoding{latin1}
\name{twinSIR_methods}
\alias{print.twinSIR}
\alias{summary.twinSIR}
\alias{AIC.twinSIR}
\alias{extractAIC.twinSIR}
\alias{vcov.twinSIR}
\alias{logLik.twinSIR}
\alias{print.summary.twinSIR}

\title{
Print, Summary and Extraction Methods for \code{"twinSIR"} Objects
}

\description{
Besides \code{print} and \code{summary} methods there are also some standard
extraction methods defined for objects of class \code{"twinSIR"}:
\code{vcov}, \code{logLik} and especially \code{AIC} and
\code{extractAIC}, which extract Akaike's Information Criterion.  Note that
special care is needed, when fitting models with parameter constraints such as
the epidemic effects \eqn{\alpha} in \code{twinSIR} models.  Parameter
constraints reduce the average increase in the maximized loglikelihood - thus
the penalty for constrained parameters should be smaller than the factor 2 used
in the ordinary definition of AIC.  To this end, these two methods offer the
calculation of the so-called one-sided AIC (OSAIC).
}

\usage{
\method{print}{twinSIR}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{twinSIR}(object,
        correlation = FALSE, symbolic.cor = FALSE, ...)

\method{AIC}{twinSIR}(object, ..., k = 2, one.sided = NULL, nsim = 1e3)
\method{extractAIC}{twinSIR}(fit, scale = 0, k = 2, one.sided = NULL,
           nsim = 1e3, ...)

\method{vcov}{twinSIR}(object, ...)
\method{logLik}{twinSIR}(object, ...)

\method{print}{summary.twinSIR}(x,
      digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), ...)
}

\arguments{
  \item{x, object, fit}{an object of class \code{"twinSIR"}.\cr
    For the \code{print} method of the \code{summary} method, an object of
    class \code{"summary.twinSIR"}.}
  \item{digits}{
    integer, used for number formatting with \code{signif()}.  Minimum number of
    significant digits to be printed in values.
  }
  \item{correlation}{
    logical. if \code{TRUE}, the correlation matrix of the estimated parameters
    is returned and printed.
  }
  \item{symbolic.cor}{
    logical. If \code{TRUE}, print the correlations in a symbolic form (see
    \code{symnum}) rather than as numbers.
  }
  \item{\dots}{
    For the \code{summary} method: arguments passed to
    \code{\link{extractAIC.twinSIR}}.\cr
    For the \code{AIC} method, optionally more fitted model objects.\cr
    For the \code{print}, \code{extractAIC}, \code{vcov} and
    \code{logLik} methods: unused (argument of the generic).
  }
  \item{k}{
    numeric specifying the "weight" of the \emph{penalty} to be used;
    in an unconstrained fit \code{k = 2} is the classical AIC.
  }
  \item{one.sided}{
    logical or \code{NULL} (the default).  Determines if the one-sided AIC
    should be calculated instead of using the classical penalty \code{k*edf}. 
    The default value \code{NULL} chooses classical AIC in the case of an
    unconstrained fit and one-sided AIC in the case of constraints. The type of
    the fit can be seen in \code{object$method} (or \code{fit$method}
    respectively), where \code{"L-BFGS"} means constrained optimization.
  }
  \item{nsim}{
    when there are more than two epidemic covariates in the fit, the
    weights in the OSAIC formula have to be determined by simulation. 
    Default is to use 1000 samples. Note that package \pkg{quadprog} is
    additionally required in this case.
  }
  \item{scale}{unused (argument of the generic).}
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance
      stars} are printed for each coefficient.}
}

\details{
  The \code{print} and \code{summary} methods allow the compact or comprehensive
  representation of the fitting results, respectively.  The former only prints
  the original function call, the estimated coefficients and the maximum
  log-likelihood value.  The latter prints the whole coefficient matrix with
  standard errors, z- and p-values (see \code{\link{printCoefmat}}), and 
  additionally the number of infections per log-baseline \code{interval}, 
  the (one-sided) AIC and the number of log-likelihood evaluations.  They both
  append a big \dQuote{WARNING}, if the optimization algorithm did not converge.

  The estimated coefficients may be extracted by using the default
  \code{coef}-method from package \pkg{stats}.
    
  The two AIC functions differ only in that \code{AIC} can take more than one
  fitted model object and that \code{extractAIC} always returns the number of
  parameters in the model (\code{AIC} only does with more than one fitted model
  object).
  
  Concerning the choice of one-sided AIC: parameter constraints -- such as the
  non-negative constraints for the epidemic effects alpha in \code{twinSIR}
  models -- reduce the average increase in the maximized loglikelihood.  Thus,
  the penalty for constrained parameters should be smaller than the factor 2
  used in the ordinary definition of AIC.  One-sided AIC (OSAIC) suggested by
  Hughes and King (2003) is such a proposal when \eqn{p} out of \eqn{k = p + q}
  parameters have non-negative constraints:
  \deqn{OSAIC = -2 l(\theta, \tau) + 2 \sum_{g=0}^p w(p,g) (k-p+g)}{%
        OSAIC = -2 l(theta, tau) + 2 sum_{g=0}^p w(p,g) (k-p+g)}
  where \eqn{w(p,g)} are \eqn{p}-specific weights.  For more details see
  Section 5.2 in \enc{Höhle}{Hoehle} (2009).
}

\value{
  The \code{print} methods return their first argument, invisibly, as
  they always should. The \code{vcov} and \code{logLik}
  methods return the estimated variance-covariance
  matrix of the parameters (here, the inverse of the estimate of the
  expected Fisher information matrix), and the maximum log-likelihood
  value of the model, respectively.
  The \code{summary} method returns a list containing some summary
  statistics of the fitted model, which is nicely printed by the
  corresponding \code{print} method.
  For the \code{\link{AIC}} and \code{\link{extractAIC}} methods, see
  the documentation of the corresponding generic functions.
}

\references{
  Hughes A, King M (2003)
  Model selection using AIC in the presence of one-sided information.
  \emph{Journal of Statistical Planning and Inference} \strong{115},
  pp. 397--411.

  \enc{Höhle}{Hoehle}, M. (2009),
  Additive-Multiplicative Regression Models for Spatio-Temporal
  Epidemics, Biometrical Journal, 51(6):961-978.
}

\author{
Michael \enc{Höhle}{Hoehle} and Sebastian Meyer
}

\examples{
data("hagelloch")

# a simplistic twinSIR model
fit <- twinSIR(~ household + cox(AGE), data = hagelloch)

coef(fit)
vcov(fit)
logLik(fit)

summary(fit, correlation = TRUE, symbolic.cor = TRUE)

# AIC or OSAIC
AIC(fit)
AIC(fit, one.sided = FALSE)
extractAIC(fit)
extractAIC(fit, one.sided = FALSE)

# comparing models via AIC
fit2 <- update(fit, nIntervals = 2)
AIC(fit, fit2)   # the 2nd column should be named "OSAIC" here
}

\keyword{methods}
\keyword{print}
\keyword{htest}
