\name{formatPval}
\alias{formatPval}

\title{
Pretty p-Value Formatting
}

\description{
  Just \acronym{yapf} -- yet another p-value formatter...
  
  It is implemented as
  \code{sapply(pv,
    function(p) format.pval(p, digits = if (p<10*eps) 1 else 2, eps = eps)}.
}

\usage{
formatPval(pv, eps = 1e-4)
}

\arguments{
  \item{pv}{a numeric vector (of p-values).}
  \item{eps}{a numerical tolerance, see \code{\link{format.pval}}.}
}

\value{
  The character vector of formatted p-values.
}

\examples{
formatPval(c(0.13567, 0.0432, 0.000546, 1e-8))
}

\keyword{print}
