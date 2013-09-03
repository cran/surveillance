\name{qlomax}
\alias{qlomax}
\title{
Quantile Function of the Lomax Distribution
}

\description{
Quantile function of the Lomax distribution with positive scale parameter
\code{scale} (often denoted as \eqn{\sigma}{sigma}) and positive shape
parameter \code{shape} (often denoted as \eqn{\alpha}{alpha}).
This implementation does not include any checks, but only the raw
formula \code{scale * ((1-p)^(-1/shape) - 1)}.
Another implementation can be found as \code{\link[VGAM]{qlomax}} in
the package \pkg{VGAM}.
}

\usage{
qlomax(p, scale, shape)
}

\arguments{
  \item{p}{
    vector of probabilities.
  }
  \item{scale}{
    positive scale parameter.
  }
  \item{shape}{
    positive shape parameter.
  }
}

\value{
  Numeric vector of quantiles corresponding to the probabilities \code{p}.
}

\author{
  Sebastian Meyer
}

\seealso{
\code{\link[VGAM]{Lomax}} in package \pkg{VGAM}.
}

\examples{
qlomax(0.99, 1, 2)
}

\keyword{distribution}
