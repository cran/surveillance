\name{all.equal}
\alias{all.equal.twinstim}
\alias{all.equal.hhh4}
\title{
Test if Two Model Fits are (Nearly) Equal
}

\description{
Two model fits are compared using standard \code{\link{all.equal}}-methods
after discarding certain elements considered irrelevant for the equality
of the fits, e.g., the runtime and the call.
}

\usage{
\method{all.equal}{twinstim}(target, current, ..., ignore = NULL)

\method{all.equal}{hhh4}(target, current, ..., ignore = NULL)
}

\arguments{
  \item{target,current}{the model fits to be compared.}
  \item{\dots}{further arguments for standard
    \code{\link{all.equal}}-methods, e.g., the numerical
    \code{tolerance}.}
  \item{ignore}{an optional character vector of elements to ignore when
    comparing the two fitted objects. The following elements are always
    ignored: \code{"runtime"} and \code{"call"}.}
}

\value{
Either \code{TRUE} or a character vector describing differences between
the \code{target} and the \code{current} model fit.
}

\author{
Sebastian Meyer
}

\keyword{utilities}
