\name{formatPval}
\alias{formatPval}

\title{
Pretty p-Value Formatting
}

\description{
  Just \acronym{yapf} -- yet another p-value formatter...
  
  It is a wrapper around \code{\link{format.pval}},
  such that by default \code{eps = 1e-4}, \code{scientific = FALSE},
  \code{digits = if (p<10*eps) 1 else 2}, and \code{nsmall = 2}.
}

\usage{
formatPval(pv, eps = 1e-4, scientific = FALSE, ...)
}

\arguments{
  \item{pv}{a numeric vector (of p-values).}
  \item{eps}{a numerical tolerance, see \code{\link{format.pval}}.}
  \item{scientific}{see \code{\link{format}}.}
  \item{\dots}{further arguments passed to \code{\link{format.pval}}
    (but \code{digits} and \code{nsmall} are hard-coded internally).}
}

\value{
  The character vector of formatted p-values.
}

\examples{
formatPval(c(0.9, 0.13567, 0.0432, 0.000546, 1e-8))
}

\keyword{print}
