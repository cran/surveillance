\name{clapply}
\alias{clapply}

\title{
Conditional \code{lapply}
}

\description{
Use \code{\link{lapply}} if the input is a list and otherwise apply the
function directly to the input \emph{and} wrap the result in a list.
The function is implemented as
\preformatted{
    if (is.list(X)) lapply(X, FUN, ...) else list(FUN(X, ...))
}
}

\usage{
clapply(X, FUN, ...)
}

\arguments{
  \item{X}{a list or a single \code{R} object on which to apply \code{FUN}.}
  \item{FUN}{the function to be applied to (each element of) \code{X}.}
  \item{\dots}{optional arguments to \code{FUN}.}
}

\value{
a list (of length 1 if \code{X} is not a list).
}

\keyword{iteration}
\keyword{list}
