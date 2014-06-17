\name{twinstim_tiaf}
\alias{tiaf}
\title{
  Temporal Interaction Function Objects
}

\description{
  A temporal interaction function for use in \code{\link{twinstim}}
  can be constructed via the \code{tiaf} function.
  It checks the supplied function elements, assigns defaults for
  missing arguments, and returns all checked arguments in a list.
  However, for standard applications it is much easier to use one of the
  pre-defined temporal interaction functions, e.g.,
  \code{\link{tiaf.exponential}}.
}

\usage{
tiaf(g, G, deriv, Deriv, npars, validpars = NULL)
}

\arguments{
  \item{g}{the temporal interaction function. It must accept
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
  \code{\link{tiaf.exponential}} for a pre-defined temporal interaction
  function, and \code{\link{siaf}} for the spatial interaction function.
}

\keyword{utilities}
