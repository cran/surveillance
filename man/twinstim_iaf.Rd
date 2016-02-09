\encoding{latin1}
\name{twinstim_iaf}
\alias{siaf.constant}
\alias{siaf.step}
\alias{siaf.gaussian}
\alias{siaf.powerlaw}
\alias{siaf.powerlawL}
\alias{siaf.student}
\alias{tiaf.constant}
\alias{tiaf.step}
\alias{tiaf.exponential}

\title{
  Temporal and Spatial Interaction Functions for \code{twinstim}
}

\description{
  A \code{twinstim} model as described in Meyer et al. (2012) requires
  the specification of the spatial and temporal interaction functions
  (\eqn{f} and \eqn{g}, respectively), 
  i.e. how infectivity decays with increasing spatial and temporal
  distance from the source of infection.
  It is of course possible to define own functions (see
  \code{\link{siaf}} and \code{\link{tiaf}}, respectively), but the
  package already predefines some useful dispersal kernels returned by
  the constructor functions documented here.
  See Meyer and Held (2014) for various spatial interaction functions,
  and Meyer et al. (2016, Section 3) for an illustration of
  the implementation.
}

\usage{
# predefined spatial interaction functions
siaf.constant()
siaf.step(knots, maxRange = Inf, nTypes = 1, validpars = NULL)
siaf.gaussian(nTypes = 1, logsd = TRUE, density = FALSE,
              F.adaptive = TRUE, effRangeMult = 6, validpars = NULL)
siaf.powerlaw(nTypes = 1, validpars = NULL)
siaf.powerlawL(nTypes = 1, validpars = NULL)
siaf.student(nTypes = 1, validpars = NULL)

# predefined temporal interaction functions
tiaf.constant()
tiaf.step(knots, maxRange = Inf, nTypes = 1, validpars = NULL)
tiaf.exponential(nTypes = 1, validpars = NULL)
}

\arguments{
  \item{knots}{numeric vector of distances at which the step function
    switches to a new height. The length of this vector determines the
    number of parameters to estimate. For identifiability, the step
    function has height 1 in the first interval \eqn{[0,knots_1)}. Note
    that the implementation is right-continuous, i.e., intervals are
    \eqn{[a,b)}.}
  \item{maxRange}{a scalar larger than any of \code{knots}.
    Per default (\code{maxRange=Inf}), the step function
    never drops to 0 but keeps the last height for any distance larger
    than the last knot. However, this might not work in some cases,
    where the last parameter value would become very small and lead to
    numerical problems. It is then possible to truncate
    interaction at a distance \code{maxRange} (just like what the
    variables \code{eps.s} and \code{eps.t} do in the
    \code{"\link{epidataCS}"} object).}
  \item{nTypes}{
    determines the number of parameters ((log-)scales or (log-)shapes)
    of the kernels. In a multitype epidemic, the different types may
    share the same spatial interaction function, in which case
    \code{nTypes=1}. Otherwise \code{nTypes} should equal the number of
    event types of the epidemic, in which case every type has its own
    (log-)scale or (log-)shape, respectively.\cr
    Currently, \code{nTypes > 1} is only implemented for
    \code{siaf.gaussian}, \code{tiaf.step}, and \code{tiaf.exponential}.
  }
  \item{logsd}{
    logical indicating if the kernel should be parametrized
    with the log-standard deviation as the parameter in question to
    enforce positivity. This is the recommended default and avoids
    constrained optimisation (L-BFGS-B) or the use of
    \code{validpars}.\cr
    The power-law kernels always use the log-scale for their scale and
    shape parameters.
  }
  \item{density}{
    logical indicating if the density or just its kernel should be used.
    If \code{density=TRUE}, \code{siaf.gaussian} uses the density of the
    bivariate, isotropic normal distribution as the spatial interaction
    function. Otherwise (default), only the kernel of the bivariate
    normal density is used.
  }
  \item{F.adaptive}{
    If \code{F.adaptive = TRUE}, then an adaptive bandwidth of
    \code{adapt*sd} will be used in the midpoint-cubature
    (\code{\link[polyCub]{polyCub.midpoint}} in package \pkg{polyCub})
    of the Gaussian interaction
    kernel, where \code{adapt} is an extra parameter of the returned
    \code{siaf$F} function and defaults to 0.1. It can be customized
    either by the \code{control.siaf$F} argument list of
    \code{twinstim}, or by a numeric specification of \code{F.adaptive}
    in the constructing call, e.g., \code{F.adaptive = 0.05} to achieve
    higher accuracy.\cr
    Otherwise, if \code{F.adaptive = FALSE}, a general cubature
    method (\code{\link[polyCub]{polyCub}}) is returned as \code{siaf$F}
    component, where the \code{method} and accuracy (\code{eps},
    \code{dimyx}, or \code{nGQ}, depending on the method) should then be
    specified in \code{twinstim}'s \code{control.siaf$F} argument.
  }
  \item{effRangeMult}{
    determines the effective range for numerical integration
    in terms of multiples of the standard deviation \eqn{\sigma} of the
    Gaussian kernel, i.e. with \code{effRangeMult=6}
    the \eqn{6 \sigma} region around the event is considered as
    the relevant integration domain instead
    of the whole observation region \code{W}.
    Setting \code{effRangeMult=NULL} will disable
    the integral approximation with an effective integration range.
  }
  \item{validpars}{
    function taking one argument, the parameter vector, indicating if it
    is valid (see also \code{\link{siaf}}).
    If \code{logsd=FALSE} and one prefers not to use
    \code{method="L-BFGS-B"} for fitting the \code{twinstim}, then
    \code{validpars} could be set to \code{function (pars) pars > 0}.
  }
}

\details{
  Evaluation of \code{twinstim}'s likelihood involves cubature of the
  spatial interaction function over polygonal domains. Various
  approaches have been compared by Meyer (2010, Section 3.2) and a new
  efficient method, which takes advantage of the assumed isotropy, has
  been proposed by Meyer and Held (2014, Supplement B, Section 2) for
  evaluation of the power-law kernels.
  These cubature methods are available in the dedicated \R package
  \pkg{polyCub} and used by the kernels implemented in \pkg{surveillance}.
  
  The readily available spatial interaction functions are defined as
  follows:
  \describe{
    \item{\code{siaf.constant}:}{
      \eqn{f(s) = 1}
    }
    \item{\code{siaf.step}:}{
      \eqn{f(s) = \sum_{k=0}^K \exp(\alpha_k) I_k(||s||)},\cr
      where \eqn{\alpha_0 = 0}, and \eqn{\alpha_1, \dots, \alpha_K} are
      the parameters (heights) to estimate. \eqn{I_k(||s||)} indicates
      if distance \eqn{||s||} belongs to the \eqn{k}th interval
      according to \code{c(0,knots,maxRange)}, where \eqn{k=0} indicates
      the interval \code{c(0,knots[1])}.\cr
      Note that \code{siaf.step} makes use of the \pkg{memoise} package
      if it is available -- and that is highly recommended to speed up
      calculations. Specifically, the areas of the intersection of a
      polygonal domain (influence region) with the \dQuote{rings} of the
      two-dimensional step function will be cached such that they are
      only calculated once for every \code{polydomain} (in the first
      iteration of the \code{twinstim} optimization). They are used in
      the integration components \code{F} and \code{Deriv}.
      See Meyer and Held (2014) for a use case and further details.
    }
    \item{\code{siaf.gaussian}:}{
      \eqn{f(s|\kappa) = \exp(-||s||/2/\sigma_\kappa^2)}\cr
      If \code{nTypes=1} (single-type epidemic or type-invariant
      \code{siaf} in multi-type epidemic), then
      \eqn{\sigma_\kappa = \sigma} for all types \eqn{\kappa}.
      If \code{density=TRUE}, then the kernel formula above is
      additionally divided by \eqn{2 \pi \sigma_\kappa^2}, yielding the
      density of the bivariate, isotropic Gaussian distribution with
      zero mean and covariance matrix \eqn{\sigma_\kappa^2 I_2}.
    }
    \item{\code{siaf.powerlaw}:}{
      \eqn{f(s) = (||s|| + \sigma)^{-d}},\cr
      which is the kernel of the Lomax density, i.e. without any
      proportionality constants. The parameters are optimized on the
      log-scale to ensure positivity, i.e.
      \eqn{\sigma = \exp(\tilde{\sigma})} and \eqn{d = \exp(\tilde{d})},
      where \eqn{(\tilde{\sigma}, \tilde{d})} is the parameter vector.
    }
    \item{\code{siaf.powerlawL}:}{
      \eqn{f(s) = (||s||/\sigma)^{-d}}, for \eqn{||s|| \ge \sigma}, and
      \eqn{f(s) = 1} otherwise,\cr
      which is a \emph{L}agged power-law kernel featuring uniform
      short-range dispersal (up to distance \eqn{\sigma}) and a
      power-law decay (Pareto-style) from distance \eqn{\sigma} onwards.
      The parameters are optimized on the log-scale to ensure positivity, i.e.
      \eqn{\sigma = \exp(\tilde{\sigma})} and \eqn{d = \exp(\tilde{d})},
      where \eqn{(\tilde{\sigma}, \tilde{d})} is the parameter vector.
      However, there is a caveat associated with this kernel: Its
      derivative wrt \eqn{\tilde{\sigma}} is mathematically undefined at
      the threshold \eqn{||s||=\sigma}. This local non-differentiability
      makes \code{twinstim}'s likelihood maximization sensitive wrt
      parameter start values, and is likely to cause false convergence
      warnings by \code{\link{nlminb}}. Possible work-arounds are to use
      the slow and robust \code{method="Nelder-Mead"}, or to just ignore
      the warning and verify the result by sets of different start values.
    }
    \item{\code{siaf.student}:}{
      \eqn{f(s) = (||s||^2 + \sigma^2)^{-d}},\cr
      which is a reparametrized \eqn{t}-kernel.
      For \eqn{d=1}, this is the kernel of the Cauchy density with scale
      \code{sigma}. In Geostatistics, a correlation function of this
      kind is known as the Cauchy model.\cr
      The parameters are optimized on the log-scale to ensure
      positivity, i.e. \eqn{\sigma = \exp(\tilde{\sigma})} and
      \eqn{d = \exp(\tilde{d})}, where \eqn{(\tilde{\sigma}, \tilde{d})}
      is the parameter vector.
    }
  }
  
  The predefined temporal interaction functions are defined as follows:
  \describe{
    \item{\code{tiaf.constant}:}{
      \eqn{g(t) = 1}
    }
    \item{\code{tiaf.step}:}{
      \eqn{g(t) = \sum_{k=0}^K \exp(\alpha_k) I_k(t)},\cr
      where \eqn{\alpha_0 = 0}, and \eqn{\alpha_1, \dots, \alpha_K} are
      the parameters (heights) to estimate. \eqn{I_k(t)} indicates
      if \eqn{t} belongs to the \eqn{k}th interval
      according to \code{c(0,knots,maxRange)}, where \eqn{k=0} indicates
      the interval \code{c(0,knots[1])}.
    }
    \item{\code{tiaf.exponential}:}{
      \eqn{g(t|\kappa) = \exp(-\alpha_\kappa t)},\cr
      which is the kernel of the exponential distribution.
      If \code{nTypes=1} (single-type epidemic or type-invariant
      \code{tiaf} in multi-type epidemic), then
      \eqn{\alpha_\kappa = \alpha} for all types \eqn{\kappa}.
    }
  }
}

\value{
  The specification of an interaction function, which is a list.
  See \code{\link{siaf}} and \code{\link{tiaf}}, respectively, for a
  description of its components.  
}

\references{
  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{ä}{ae}t
  M\enc{ü}{ue}nchen.\cr
  Available as \url{http://epub.ub.uni-muenchen.de/11703/}

  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}

  Meyer, S. and Held, L. (2014):
  Power-law models for infectious disease spread.
  \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.\cr
  DOI-Link: \url{http://dx.doi.org/10.1214/14-AOAS743}

  Meyer, S., Held, L. and H\enc{ö}{oe}hle, M. (2016):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}. \emph{Journal of Statistical Software}. In press.\cr
  Preprint available at \url{http://arxiv.org/abs/1411.0416}
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{twinstim}}, \code{\link{siaf}}, \code{\link{tiaf}},
  and package \pkg{polyCub} for the involved cubature methods.
}

\examples{
# constant temporal dispersal
tiaf.constant()
# step function kernel
tiaf.step(c(3,7), maxRange=14, nTypes=2)
# exponential decay specification
tiaf.exponential()

# Type-dependent Gaussian spatial interaction function using an adaptive
# two-dimensional midpoint-rule to integrate it over polygonal domains
siaf.gaussian(2, F.adaptive=TRUE)

# Type-independent power-law kernel
siaf.powerlaw()

# "lagged" power-law
siaf.powerlawL()

# (reparametrized) t-kernel
siaf.student()

# step function kernel
siaf.step(c(10,20,50), maxRange=100)
}

\keyword{models}
\keyword{utilities}
