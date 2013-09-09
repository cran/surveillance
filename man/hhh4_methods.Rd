\encoding{latin1}
\name{hhh4_methods}
\alias{print.ah4}
\alias{summary.ah4}
\alias{logLik.ah4}
\alias{AIC.ah4}
\alias{coef.ah4}
\alias{ranef}
\alias{fixef}
\alias{fixef.ah4}
\alias{ranef.ah4}
\alias{confint.ah4}
% FIXME: functions without documentation atm
%% \alias{print.summary.ah4}
%% \alias{plot.ah4}
%% \alias{update.ah4}

\title{
Print, Summary and Extraction Methods for \code{"ah4"} Objects
}

\description{
Besides \code{print} and \code{summary} methods there are also some standard
extraction methods defined for objects of class \code{"ah4"} resulting
from a call to \code{\link{hhh4}}.
}

\usage{
\method{print}{ah4}(x, digits = max(3, getOption("digits") - 3),
       reparamPsi = TRUE, ...)
\method{summary}{ah4}(object, ...)

\method{coef}{ah4}(object, se = FALSE, reparamPsi = TRUE, 
       idx2Exp = NULL, amplitudeShift = FALSE, ...)
\method{fixef}{ah4}(object, ...)
\method{ranef}{ah4}(object, tomatrix = FALSE, ...)

\method{AIC}{ah4}(object, ..., k = 2)
\method{logLik}{ah4}(object, ...)

\method{confint}{ah4}(object, parm, level = 0.95, 
       reparamPsi = TRUE, idx2Exp = NULL, amplitudeShift = FALSE, ...)
}

\arguments{
  \item{x, object}{an object of class \code{"ah4"}.}
  \item{digits}{the number of significant digits to use when printing
  }
  \item{reparamPsi}{
    logical. If \code{TRUE}, the overdispersion parameter from the 
    negative binomial distribution is transformed from internal (log-)scale
    to a user scale (where zero corresponds to a Poisson distribution). 
  }
  \item{\dots}{
    For the \code{summary}, \code{fixef} and \code{ranef} methods: 
    arguments passed to \code{coef}.\cr
    For the remaining methods: unused (argument of the generic).
  }
  \item{se}{logical switch indicating if standard errors are required}
  \item{idx2Exp}{vector with integers idicating the parameters
    which should be returned on exp-scale.}
  \item{amplitudeShift}{logical switch indicating whether the parameters
   for sine/cosine terms modelling seasonal patterns 
   (see \code{\link{addSeason2formula}}) should be transformed
   to an amplitude/shift formulation.}
  \item{tomatrix}{logical. If \code{FALSE} (default), the vector of
    all random effects is returned (as used internally). However, for
    random intercepts of \code{type="car"}, the number of parameters is
    one less than the number of regions and the individual parameters are
    not obviously linked to specific regions. Setting \code{tomatrix} to
    \code{TRUE} returns a more useful representation of random effects
    in a matrix with as many rows as there are regions and as many
    columns as there are random effects. Here, any CAR-effects are
    transformed to region-specific effects.}
  \item{k}{
    numeric specifying the \emph{penalty} per parameter to be used;
    \code{k = 2} (default) is the classical AIC.
  }
  \item{parm}{a vector of numbers or names, specifying which parameters 
    are to be given confidence intervals. If missing, all parameters
    are considered.}
  \item{level}{the confidence level required.}
}

\details{
  The function \code{coef} returns all estimated (regression)
  parameters from a model as propsed by Paul and Held (2011), see 
  \code{\link{hhh4}}. 
  If the model includes random effects, those can be extracted with
  \code{ranef}, whereas \code{fixef} returns the fixed
  parameters.

  The function \code{AIC} returns the value of the AIC criterion
  only for models without random effects; otherwise it is not well
  defined and \code{NA} is returned.

  The function \code{logLik} returns an object of class \code{"logLik"}
  in the case of a model without random effects, and the value of the 
  penalized log-likelihood at the parameter estimates otherwise.

  The function \code{confint} returns Wald-type confidence intervals 
  (assuming asymptotic normality) 
  for one or more parameters in the fitted model.
}

\value{
  See also the documentation for the generic functions \code{\link{coef}},
  \code{\link{AIC}}, \code{\link{logLik}} and \code{\link{confint}}.
}

\author{
  Michaela Paul
}

\references{
  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136.
}

%\examples{
%}

\keyword{methods}
\keyword{print}
