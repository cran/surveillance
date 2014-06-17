\encoding{latin1}
\name{hhh4_methods}
\alias{print.hhh4}
\alias{summary.hhh4}
\alias{nobs.hhh4}
\alias{logLik.hhh4}
\alias{coef.hhh4}
\alias{vcov.hhh4}
\alias{ranef}
\alias{fixef}
\alias{fixef.hhh4}
\alias{ranef.hhh4}
\alias{confint.hhh4}
% FIXME: functions without documentation atm
%% \alias{print.summary.hhh4}
%% \alias{update.hhh4}
%% \alias{terms.hhh4}
%% \alias{residuals.hhh4}

\title{
Print, Summary and other Standard Methods for \code{"hhh4"} Objects
}

\description{
Besides \code{print} and \code{summary} methods there are also some standard
extraction methods defined for objects of class \code{"hhh4"} resulting
from a call to \code{\link{hhh4}}.
}

\usage{
\method{print}{hhh4}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{hhh4}(object, maxEV = FALSE, ...)

\method{coef}{hhh4}(object, se = FALSE, reparamPsi = TRUE, 
     idx2Exp = NULL, amplitudeShift = FALSE, ...)
\method{fixef}{hhh4}(object, ...)
\method{ranef}{hhh4}(object, tomatrix = FALSE, ...)

\method{nobs}{hhh4}(object, ...)
\method{logLik}{hhh4}(object, ...)

\method{vcov}{hhh4}(object, reparamPsi = TRUE, 
     idx2Exp = NULL, amplitudeShift = FALSE, ...)
\method{confint}{hhh4}(object, parm, level = 0.95, 
        reparamPsi = TRUE, idx2Exp = NULL, amplitudeShift = FALSE, ...)
}

\arguments{
  \item{x, object}{an object of class \code{"hhh4"}.}
  \item{digits}{the number of significant digits to use when printing
  }
  \item{maxEV}{logical indicating if the summary should contain the
    (range of the) dominant eigenvalue as a measure of the importance of
    the epidemic components. By default, the value is not calculated as
    this may take some seconds depending on the number of time points
    and units in \code{object$stsObj}.}
  \item{\dots}{
    For the \code{print}, \code{summary}, \code{fixef} and \code{ranef}
    methods: arguments passed to \code{coef}.\cr
    For the remaining methods: unused (argument of the generic).
  }
  \item{reparamPsi}{
    logical. If \code{TRUE} (default), the overdispersion parameter from the 
    negative binomial distribution is transformed from internal scale (-log)
    to standard scale, where zero corresponds to a Poisson distribution. 
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
  \item{parm}{a vector of numbers or names, specifying which parameters 
    are to be given confidence intervals. If missing, all parameters
    are considered.}
  \item{level}{the confidence level required.}
}

\value{
  The \code{\link{coef}}-method returns all estimated (regression)
  parameters from a \code{\link{hhh4}} model as proposed by Paul and
  Held (2011). 
  If the model includes random effects, those can be extracted with
  \code{ranef}, whereas \code{fixef} returns the fixed parameters.

  The \code{\link{nobs}}-method returns the number of observations used
  for model fitting.
  The \code{\link{logLik}}-method returns an object of class
  \code{"logLik"} with \code{"df"} and \code{"nobs"} attributes.
  For a random effects model, the value of the \emph{penalized}
  log-likelihood at the MLE is returned, but degrees of freedom are
  not available (\code{NA_real_}).
  As a consequence, \code{\link{AIC}} and \code{\link{BIC}} are only
  well defined for models without random effects;
  otherwise these functions return \code{NA_real_}.

  The \code{\link{vcov}}-method returns the estimated
  variance-covariance matrix of the \emph{regression} parameters.
  The estimated variance-covariance matrix of random effects is
  available as \code{object$Sigma}.
  The \code{\link{confint}}-method returns Wald-type confidence
  intervals (assuming asymptotic normality).
}

\seealso{
  the \code{\link[=plot.hhh4]{plot}}-method for fitted \code{hhh4} models.
}

\author{
  Michaela Paul and Sebastian Meyer
}

\references{
  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136.
}

\keyword{methods}
\keyword{print}
