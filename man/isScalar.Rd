\name{isScalar}
\alias{isScalar}

\title{
Checks if the Argument is Scalar
}

\description{
The simple helper function \code{isScalar} just checks if its argument is
a scalar, i.e. a numeric vector of length 1.
It is implemented as \code{length(x) == 1L && is.vector(x, mode = "numeric")}.
}

\usage{
isScalar(x)
}

\arguments{
  \item{x}{an \code{R} object.}
}

\value{
A length-one logical vector.
}

%% \examples{
%% isScalar(TRUE)   # FALSE
%% isScalar(1:10)   # FALSE
%% isScalar(pi)     # TRUE
%% }

\keyword{internal}
