\encoding{latin1}
\name{twinstim_iaf}
\alias{siaf.constant}
\alias{siaf.step}
\alias{siaf.gaussian}
\alias{siaf.exponential}
\alias{siaf.powerlaw}
\alias{siaf.powerlaw1}
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
  Own such functions can be specified (see
  \code{\link{siaf}} and \code{\link{tiaf}}, respectively), but the
  package already predefines some common dispersal kernels returned by
  the constructor functions documented here.
  See Meyer and Held (2014) for various spatial interaction functions,
  and Meyer et al. (2017, Section 3, available as \code{vignette("twinstim")})
  for an illustration of the implementation.
}

\usage{
# predefined spatial interaction functions
siaf.constant()
siaf.step(knots, maxRange = Inf, nTypes = 1, validpars = NULL)
siaf.gaussian(nTypes = 1, logsd = TRUE, density = FALSE,
              F.adaptive = FALSE, F.method = "iso",
              effRangeMult = 6, validpars = NULL)
siaf.exponential(nTypes = 1, validpars = NULL, engine = "C")
siaf.powerlaw(nTypes = 1, validpars = NULL, engine = "C")
siaf.powerlaw1(nTypes = 1, validpars = NULL, sigma = 1)
siaf.powerlawL(nTypes = 1, validpars = NULL, engine = "C")
siaf.student(nTypes = 1, validpars = NULL, engine = "C")

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
    \eqn{[a,b)}.\cr
    An initial choice of knots could be based on quantiles of the
    observed distances between events and their potential source events.
    For instance, an identifiable spatial step function could be
    \code{siaf.step(quantile(\link{getSourceDists}(myepi, "space"), c(1,2,4)/10))},
    where \code{myepi} is the \code{"epidataCS"} data to be modelled.}
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
    \code{siaf.gaussian(F.adaptive = TRUE)},
    \code{tiaf.step}, and \code{tiaf.exponential}.
  }
  \item{logsd,density}{
    logicals affecting the parametrization of the Gaussian kernel.
    Settings different from the defaults are deprecated.
    The default is to use only the kernel of the bivariate, isotropic
    normal distribution (\code{density=FALSE}, see Details below),
    parametrized with the log-standard deviation (\code{logsd=TRUE}) to
    avoid constrained optimisation (L-BFGS-B) or \code{validpars}.\cr
    The power-law kernels always employ the log-scale for their scale
    and shape parameters.
  }
  \item{F.adaptive,F.method}{
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
    Otherwise, if \code{F.adaptive = FALSE}, the \code{F.method}
    argument determines which \code{\link[polyCub]{polyCub}} method to
    use in \code{siaf$F}. The accuracy (controlled via, e.g.,
    \code{nGQ}, \code{rel.tol}, or \code{eps}, depending on the cubature
    method) can then be adjusted in \code{twinstim}'s
    \code{control.siaf$F} argument.
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
  \item{engine}{
    character string specifying the implementation to use.
    Prior to \pkg{surveillance} 0.14.0, the \code{intrfr} functions for
    \code{\link{polyCub.iso}} were evaluated in \R (and this
    implementation is available via \code{engine = "R"}).
    The new C-implementation, \samp{LinkingTo} the newly exported
    \code{polyCub_iso} C-implementation in \pkg{polyCub} 0.6.0,
    is considerably faster.
  }
  \item{sigma}{
    Fixed value of \eqn{\sigma} for the one-parameter power-law kernel.
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
      If \code{density=TRUE} (deprecated), then the kernel formula above is
      additionally divided by \eqn{2 \pi \sigma_\kappa^2}, yielding the
      density of the bivariate, isotropic Gaussian distribution with
      zero mean and covariance matrix \eqn{\sigma_\kappa^2 I_2}.
      The standard deviation is optimized on the log-scale
      (\code{logsd = TRUE}, not doing so is deprecated).
    }
    \item{\code{siaf.exponential}:}{
      \eqn{f(s) = exp(-||s||/sigma)}\cr
      The scale parameter \eqn{sigma} is estimated on the log-scale,
      i.e., \eqn{\sigma = \exp(\tilde{\sigma})}, and \eqn{\tilde{\sigma}}
      is the actual model parameter.
    }
    \item{\code{siaf.powerlaw}:}{
      \eqn{f(s) = (||s|| + \sigma)^{-d}}\cr
      The parameters are optimized on the log-scale to ensure positivity, i.e.,
      \eqn{\sigma = \exp(\tilde{\sigma})} and \eqn{d = \exp(\tilde{d})},
      where \eqn{(\tilde{\sigma}, \tilde{d})} is the parameter vector.
      If a power-law kernel is not identifiable for the dataset at hand,
      the exponential kernel or a lagged power law are useful alternatives.
    }
    \item{\code{siaf.powerlaw1}:}{
      \eqn{f(s) = (||s|| + 1)^{-d}},\cr
      i.e., \code{siaf.powerlaw} with fixed \eqn{\sigma = 1}.
      A different fixed value for \eqn{sigma} can be specified via the
      \code{sigma} argument of \code{siaf.powerlaw1}.
      The decay parameter \eqn{d} is estimated on the log-scale.
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
      warnings by \code{\link{nlminb}}. Possible workarounds are to use
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
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}

  Meyer, S. and Held, L. (2014):
  Power-law models for infectious disease spread.
  \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.
  \doi{10.1214/14-AOAS743}

  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
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
# exponential temporal decay
tiaf.exponential()

# Type-dependent Gaussian spatial interaction function using an adaptive
# two-dimensional midpoint-rule to integrate it over polygonal domains
siaf.gaussian(2, F.adaptive=TRUE)

# Single-type Gaussian spatial interaction function (using polyCub.iso)
siaf.gaussian()

# Exponential kernel
siaf.exponential()

# Power-law kernel
siaf.powerlaw()

# Power-law kernel with fixed sigma = 1
siaf.powerlaw1()

# "lagged" power-law
siaf.powerlawL()

# (reparametrized) t-kernel
siaf.student()

# step function kernel
siaf.step(c(10,20,50), maxRange=100)
}

\keyword{models}
\keyword{utilities}
