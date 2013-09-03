\encoding{latin1}
\name{twinSIR_profile}
\alias{profile.twinSIR}

\title{
Profile Likelihood Computation and Confidence Intervals
}

\description{
  Function to compute estimated and profile likelihood based confidence
  intervals. Computations might be cumbersome!
}

\usage{
\method{profile}{twinSIR}(fitted, profile, alpha = 0.05,
        control = list(fnscale = -1, factr = 10, maxit = 100), ...)
}

\arguments{
  \item{fitted}{
    an object of class \code{"twinSIR"}.
  }
  \item{profile}{
    a list with elements being numeric vectors of length 4.  These vectors must
    have the form \code{c(index, lower, upper, gridsize)}.
    \describe{
      \item{\code{index}:}{
        index of the parameter to be profiled in the vector \code{coef(fitted)}.
      }
      \item{\code{lower, upper}:}{
        lower/upper limit of the grid on which the profile log-likelihood is
        evaluated. Can also be \code{NA} in which case \code{lower/upper} equals 
        the lower/upper bound of the respective 0.3 \% Wald confidence interval
        (+-3*se).
      }
      \item{\code{gridsize}:}{
        grid size of the equally spaced grid between lower and upper.  Can also
        be 0 in which case the profile log-likelihood for this parameter is not 
        evaluated on a grid.
      }
    }
  }
  \item{alpha}{
    \eqn{(1-\alpha)\%}{(1-alpha)\%} profile likelihood based confidence
    intervals are computed.  If alpha <= 0, then no confidence intervals are
    computed.
  }
  \item{control}{
    control object to use in \code{\link{optim}} for the profile log-likelihood
    computations.
  }
  \item{\dots}{
    unused (argument of the generic).
  }
}

\value{
list with profile log-likelihood evaluations on the grid and highest likelihood
and wald confidence intervals.  The argument \code{profile} is also returned.
}

\author{
Michael \enc{Höhle}{Hoehle}
}

\examples{
data("foofit")
# the following call takes a while
\dontrun{
prof <- profile(foofit, list(c(1,0,0.05,5), c(3,NA,NA,0), c(4, NA, NA, 10)))
prof
}
}

\keyword{htest}
\keyword{methods}
\keyword{optimize}
\keyword{dplot}
