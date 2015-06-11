\name{coeflist}
\alias{coeflist}
\alias{coeflist.default}
\title{
  List Coefficients by Model Component
}

\description{
  S3-generic function to use with models which contain several groups of
  coefficients in their coefficient vector. The \code{coeflist} methods
  are intended to list the coefficients by group. The default method
  simply \code{\link{split}}s the coefficient vector given the number of
  coefficients by group.
}

\usage{
coeflist(x, ...)

\method{coeflist}{default}(x, npars, ...)
}

\arguments{
  \item{x}{
    a model with groups of coefficients or, for the default method, a
    vector of coefficients.
  }
  \item{npars}{
    a named vector specifying the number of coefficients per group.
  }
  \item{\dots}{
    potential further arguments (currently ignored).
  }
}

\value{
  a list of coefficients
}

\author{
  Sebastian Meyer
}

\examples{
## the default method just 'split's the coefficient vector
coefs <- c(a = 1, b = 3, dispersion = 0.5)
npars <- c(regression = 2, variance = 1)
coeflist(coefs, npars)
}

\keyword{models}
\keyword{utilities}
