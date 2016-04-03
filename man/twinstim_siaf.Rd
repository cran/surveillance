\name{twinstim_siaf}
\alias{siaf}
\title{
  Spatial Interaction Function Objects
}

\description{
  A spatial interaction function for use in \code{\link{twinstim}}
  can be constructed via the \code{siaf} function.
  It checks the supplied function elements, assigns defaults for
  missing arguments, and returns all checked arguments in a list.
  However, for standard applications it is much easier to use one of the
  pre-defined spatial interaction functions, e.g.,
  \code{\link{siaf.gaussian}}.
}

\usage{
siaf(f, F, Fcircle, effRange, deriv, Deriv, simulate, npars,
     validpars = NULL)
}

\arguments{
  \item{f}{the spatial interaction function. It must accept
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
    the \code{control.siaf$F} argument list of \code{twinstim}. If the \code{F}
    function is missing, a general default (\code{\link[polyCub]{polyCub}}) 
    will be used, with extra arguments \code{method} (default: \code{"SV"}) 
    and corresponding accuracy parameters.}

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
    second argument) over a polygonal \code{"owin"} domain (first
    argument). The return value is thus a vector of length \code{npars}.
    The third argument is the parameter vector and the fourth argument
    is a (\emph{single}) type and must be named \code{type}.
    There may be additional arguments, which can then be specified in
    the \code{control.siaf$Deriv} argument list of \code{twinstim}. If the
    \code{Deriv} function is missing, a general default
    (\code{\link[polyCub]{polyCub}}) will be used, with extra arguments
    \code{method} (default: \code{"SV"}) and corresponding accuracy parameters.}

  \item{simulate}{optional function returning a sample drawn from the
    spatial kernel (only required for the simulation of \code{twinstim}
    models).  Its first argument is the size of the sample to
    generate, next the parameter vector, an optional single event type,
    and an optional upper bound for the radius within which to simulate
    points. The function must return a two-column \emph{matrix} of the
    sampled locations. 
    Note that the simulation method actually samples only one location
    at a time, thus it is sufficient to have a working
    \code{function(n=1, pars, type, ub)}.
  }

  \item{npars}{the number of parameters of the spatial interaction
    function \code{f} (i.e. the length of its second argument).}
  
  \item{validpars}{
    optional function taking one argument, the parameter vector, indicating if it
    is valid. This approach to specify parameter constraints is rarely
    needed, because usual box-constrained parameters can be taken into
    account by using L-BFGS-B as the optimization method in
    \code{twinstim} (with arguments \code{lower} and \code{upper}), and
    positivity constraints by using log-parametrizations.
    This component is not necessary (and ignored) if \code{npars == 0}.
  }
}

\value{
  list of checked arguments.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{siaf.gaussian}} for a pre-defined spatial interaction
  function, and \code{\link{tiaf}} for the temporal interaction function.
}

\keyword{utilities}
