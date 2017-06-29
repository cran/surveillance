\name{hhh4_internals}
\alias{meanHHH}
\alias{sizeHHH}
\alias{decompose.hhh4}
\title{
  Internal Functions Dealing with \code{hhh4} Models
}
\description{
  The functions documented here are considered \emph{internal},
  i.e., not intended to be called by the user. They are used by
  add-on packages dealing with \code{\link{hhh4}} models.
}
\usage{
meanHHH(theta, model, subset = model$subset, total.only = FALSE)
sizeHHH(theta, model, subset = model$subset)

decompose.hhh4(x, coefs = x$coefficients, ...)
}
\arguments{
  \item{theta,coefs}{numeric vector of model parameters.}
  \item{model}{the model terms as returned by the
    \code{\link{terms}}-method for \code{"hhh4"} objects.}
  \item{subset}{vector of time points for which to compute the component
    means. Defaults to the fitted time range. For \code{sizeHHH},
    \code{subset=NULL} means to return the vector of dispersion
    parameters.}
  \item{total.only}{logical. Should only the total mean (epidemic +
    endemic) be returned in a \code{length(subset)} x nUnit matrix?
    Otherwise, a list of such matrices is returned, giving the values of
    the various model components separately (as well as the total).}
  \item{x}{a fitted \code{hhh4} model.}
  \item{\dots}{unused.}
}
\details{
  \code{meanHHH} computes the components of the mean returned in
  \code{length(subset)} x nUnit matrices.
  \code{sizeHHH} computes the model dispersion in \code{\link{dnbinom}}
  (\code{mu}, \code{size}) parametrization (it returns \code{NULL} in
  the Poisson case).
  \code{decompose.hhh4} decomposes the fitted mean (extracted via
  \code{meanHHH}) in an array with dimensions \eqn{(t, i, j)}, where the
  first \eqn{j} index is \code{"endemic"}.
}
\author{
  Michaela Paul and Sebastian Meyer
}
\keyword{internal}
