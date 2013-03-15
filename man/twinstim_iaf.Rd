\encoding{latin1}
\name{twinstim_iaf}
\alias{siaf.constant}
\alias{siaf.gaussian}
\alias{siaf.lomax}
\alias{siaf.powerlaw}
\alias{tiaf.constant}
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
  It is of course possible to define own functions (which are
  implemented as list objects with required components described in the
  \dQuote{Value} section below), but the package also predefines some 
  dispersal kernels returned by the constructor functions documented
  here.
}

\usage{
# predefined spatial interaction functions
siaf.constant()
siaf.gaussian(nTypes = 1, logsd = TRUE, density = FALSE, F.adaptive = TRUE,
              effRangeMult = 6, validpars = NULL)
siaf.lomax(nTypes = 1, logpars = TRUE, density = FALSE,
           effRangeProb = NULL, validpars = NULL)
siaf.powerlaw(nTypes = 1, logpars = TRUE,
              effRangeProb = NULL, validpars = NULL)

# predefined temporal interaction functions
tiaf.constant()
tiaf.exponential(nTypes = 1)
}

\arguments{
  \item{nTypes}{
    determines the number of parameters ((log-)scales or (log-)shapes)
    of the kernels. In a multitype epidemic, the different types may
    share the same spatial interaction function, in which case
    \code{nTypes=1}. Otherwise \code{nTypes} should equal the number of
    event types of the epidemic, in which case every type has its own
    (log-)scale or (log-)shape, respectively. Currently, the constructor
    function for the Lomax and the power law kernels \code{siaf.lomax}
    and \code{siaf.powerlaw}, respectively, only allows a type-invariant
    kernel (\code{nTypes=1}).
  }
  \item{logsd, logpars}{
    logical indicating if the kernel should be parametrized
    with the log-standard deviation (or log(scale) and log(shape)) as
    the parameter(s) in question to enforce positivity. This
    is the recommended default and avoids constrained optimisation
    (L-BFGS-B) or the use of \code{validpars}. The alternative
    (\code{FALSE}) is not implemented for \code{siaf.lomax} and
    \code{siaf.powerlaw}.
  }
  \item{density}{
    logical. If \code{TRUE}, a density is used as the
    interaction function: for \code{siaf.gaussian}, the density of the
    bivariate, isotropic normal distribution, for \code{siaf.lomax},
    \eqn{f(||\bold{s}||)} is the density of the (univariate) Lomax
    distribution (\eqn{f(\bold{s})} is not a bivariate density here).
    However, in this case, the magnitude of the
    interaction function depends on the parameters of the distribution,
    which interferes with the identification of an epidemic
    intercept. Hence, the use of \code{density = TRUE} is not
    recommended.\cr
    For \code{density = FALSE}, \code{siaf.gaussian} uses the kernel of
    the bivariate normal density, and \code{siaf.lomax} is normalized to
    unity at the origin, i.e. \eqn{f(\bold{s}) / f(\bold{0})} is used.
    However, the latter still has proportionality constants in the kernel,
    which is why the use \code{siaf.lomax} is deprecated.
    The new \code{siaf.powerlaw} replaces \code{siaf.lomax} and uses the
    pure kernel of the Lomax density (see Details).
  }
  \item{F.adaptive}{
    If \code{F.adaptive = TRUE}, then an adaptive bandwidth of
    \code{adapt*sd} will be used in the midpoint-cubature
    (see \code{\link{polyCub.midpoint}}) of the Gaussian interaction
    kernel, where \code{adapt} is an extra parameter of the returned
    \code{siaf$F} function and defaults to 0.1, but it may also be
    specified in the \code{control.siaf$F} argument list of
    \code{twinstim}.\cr  
    Otherwise, if \code{F.adaptive = FALSE}, a general cubature
    method (\code{\link{polyCub}}) is returned as \code{siaf$F}
    component, where the \code{method} and accuracy (\code{eps},
    \code{dimyx}, or \code{nGQ}, depending on the method) should then be
    specified in \code{twinstim}'s \code{control.siaf$F} argument.
  }
  \item{effRangeMult}{
    determines the effective range for numerical integration
    in terms of multiples of the standard deviation \eqn{\sigma} of the
    Gaussian kernel, i.e. with \code{effRangeMult=6} numerical
    integration only considers the \eqn{6 \sigma} area around the event instead
    of the whole observation region \code{W}.
    Setting \code{effRangeMult=NULL} will disable
    the integral approximation with an effective integration range.
  }
  \item{effRangeProb}{
    numeric scalar in \eqn{(0,1)}, which determines the effective range for
    numerical integration in terms of the corresponding quantile of the
    univariate Lomax distribution. May also be \code{NULL} (default),
    which disables the integral approximation with an effective
    integration range.
  }
  \item{validpars}{
    function taking one argument, the parameter vector, indicating if it
    is valid. If \code{logsd=FALSE} or \code{logpars=FALSE} and one
    prefers not to use \code{method="L-BFGS-B"} for fitting the
    \code{twinstim}, then \code{validpars} should be set to
    \code{function (pars) pars > 0}.
  }
}

\details{
  The readily available spatial interaction functions are defined as
  follows:
  \describe{
    \item{\code{siaf.constant}:}{
      \eqn{f(s) = 1}
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
    \item{\code{siaf.lomax}:}{
      \eqn{f(s) = (\sigma / (||s|| + \sigma))^{\alpha+1}},\cr
      which is the density of the univariate Lomax distribution scaled
      to 1 at the origin, i.e. divided by \eqn{\alpha / \sigma}.
      If \code{density=TRUE}, the original Lomax density is used.
      Due to potential identifiability issues with an epidemic
      intercept, \code{siaf.lomax} is deprecated and
      \code{siaf.powerlaw} should be used instead!
    }
    \item{\code{siaf.powerlaw}:}{
      \eqn{f(s) = (||s|| + \sigma)^{-d}},\cr
      which is the pure kernel of the Lomax density, i.e. without any
      proportionality constants. The parameters are optimized on the
      log-scale to ensure positivity, i.e.
      \eqn{\sigma = \exp(\tilde{\sigma})} and \eqn{d = \exp(\tilde{d})},
      where \eqn{(\tilde{\sigma}, \tilde{d})} is the parameter vector.
    }
  }
  
  Furthermore, there are two predefined temporal interaction functions:
  \describe{
    \item{\code{tiaf.constant}:}{
      \eqn{g(t) = 1}
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
  
  For spatial kernels, this list contains the following components:
  \item{f}{spatial interaction function (\code{siaf}). It must accept
    two arguments, the first one being a (2-column) coordinate matrix, the
    second one a parameter vector. For marked \code{twinstim}, it must
    accept the type of the event (integer code) as its third argument
    (either a single type for all locations or separate types for each
    location).}
  \item{F}{function computing the integral of \eqn{f(s)} (passed as
    second argument) over a polygonal \code{"owin"} domain (first argument).
    The third and fourth argument are
    the parameter vector and the (\emph{single}) type, respectively.
    There may be additional arguments, which can then be specified in
    the \code{control.siaf$F} argument list of \code{twinstim}. If this
    function is missing, a general default (\code{polyCub}) will be
    used, with extra arguments \code{method} (default: \code{"SV"}) and
    corresponding accuracy parameters.}
  \item{Fcircle}{optional function for fast calculation of the
    (two-dimensional) integral of \eqn{f(s)} over a circle with radius
    \code{r} (first argument). Further arguments are as for \code{f}. It
    must not be vectorized (will always be called with single radius
    and a single type). If this function is specified, integration of
    the \code{siaf} over the spatial influence region of an event will
    be faster if the region is actually circular. This is the case if
    the event is located at least a distance \code{eps.s} from the border
    of the observation region \code{W}, or if the distance to the border
    is larger than the effective integration range (if specified, see
    \code{effRange} below).}
  \item{effRange}{optional function returning the \dQuote{effective}
    range of \eqn{f(s)} for the given set of parameters (the first and
    only argument) such that the circle with radius \code{effRange}
    contains the numerically essential proportion of the 
    integral mass. For the Gaussian kernel the default is
    \code{function (logsd) 6*exp(logsd)}. The return value must be a
    vector of length \code{nTypes} (effective range for each type). This
    function is only used if \code{Fcircle} is also specified.}
  \item{deriv}{optional derivative of \eqn{f(s)} \emph{with respect to
      the parameters}. It takes the same arguments as \code{f} but
    returns a matrix with as many rows as there were coordinates in the
    input and \code{npars} columns. This derivative is necessary for the
    calculation of the score function in \code{twinstim()}, which is
    advantageous for the numerical log-likelihood maximization.} 
  \item{Deriv}{function computing the integral of \code{deriv} (passed as
    second argument) over a polygonal \code{"owin"} domain (first argument).
    The third and fourth argument are
    the parameter vector and the (\emph{single}) type, respectively.
    There may be additional arguments, which can then be specified in
    the \code{control.siaf$Deriv} argument list of \code{twinstim}. If this
    function is missing, a general default (\code{polyCub}) will be
    used, with extra arguments \code{method} (default: \code{"SV"}) and
    corresponding accuracy parameters.}
  \item{simulate}{optional function returning a sample drawn from the
    spatial kernel (only required for the simulation of \code{twinstim}
    models).  Its first argument is the size of the sample to
    generate, further arguments are as for \code{f} (but the optional
    third argument specifies a single event type only).  The function must
    return a two-column \emph{matrix} of the sampled locations.
    Note that the simulation method actually samples only one location
    at a time, thus it is sufficient to have a working
    \code{function(n=1, pars, type)}. For \code{siaf.constant()}, the
    \code{simulate} component is \code{NULL}, but this will be handled
    specially by \code{simEpidataCS}.
  }
  \item{npars}{the number of parameters of the spatial interaction
    function \code{f} (i.e. the length of its second argument).}
    
  For temporal kernels, the components are:
  \item{g}{temporal interaction function (\code{tiaf}). It must accept
    two arguments, the first one being a vector of time points, the
    second one a parameter vector. For marked \code{twinstim}, it must
    accept the type of the event (integer code) as its third argument
    (either a single type for all locations or separate types for each
    location).}
  \item{G}{a primitive of \eqn{g(t)} (with respect to time). It must
    accept the same arguments as \code{g}, for instance a \emph{vector}
    of time points (not just a single one).}
  \item{deriv}{optional derivative of \eqn{g(t)} \emph{with respect to
      the parameters}. It takes the same arguments as \code{g} but
    returns a matrix with as many rows as there were time points in the
    input and \code{npars} columns. This derivative is necessary for the
    calculation of the score function in \code{twinstim()}, which is
    advantageous for the numerical log-likelihood maximization.}
  \item{Deriv}{optional primitive of \code{deriv} (with respect to
    time). It must accept the same arguments as \code{deriv}, \code{g} and
    \code{G} and returns a matrix with as many rows as there were time
    points in the input and \code{npars} columns. The integrated
    derivative is necessary for the score function in \code{twinstim}.}
  \item{npars}{the number of parameters of the temporal interaction
    function \code{g} (i.e. the length of its second argument).}
  
  In \emph{both} cases, there may be an additional component
  \code{validpars} as a way to specify parameter constraints:
  \item{validpars}{
    optional function taking one argument, the parameter vector, indicating if it
    is valid. This function is rarely needed in practice, because usual
    box-constrained parameters can be taken into account by using L-BFGS-B
    as the optimization method in \code{twinstim} (with arguments
    \code{lower} and \code{upper}), and positivity constraints by using
    log-parametrizations. This component is not necessary (and ignored) if
    \code{npars == 0}.
  }
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}

  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{ä}{ae}t
  M\enc{ü}{ue}nchen.\cr
  Available as \url{http://epub.ub.uni-muenchen.de/11703/}
}

\author{
  Sebastian Meyer
}

\seealso{
\code{\link{twinstim}}.
}

\examples{
# constant temporal dispersal
tiaf.constant()

# Type-dependent Gaussian spatial interaction function using an adaptive
# two-dimensional midpoint-rule to integrate it over polygonal domains
siaf.gaussian(2, F.adaptive=TRUE)

# Type-independent power-law kernel
siaf.powerlaw()   # siaf.lomax() is deprecated
}

\keyword{models}
\keyword{utilities}
