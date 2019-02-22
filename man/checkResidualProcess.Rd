\name{checkResidualProcess}
\alias{checkResidualProcess}
\title{
Check the residual process of a fitted \code{twinSIR} or \code{twinstim}
}

\description{
  Transform the residual process (cf. the
  \code{\link[=residuals.twinstim]{residuals}} methods for classes
  \code{"twinSIR"} and \code{"twinstim"}) such that the transformed
  residuals should be uniformly distributed if the fitted model
  well describes the true conditional intensity function. Graphically
  check this using \code{\link{ks.plot.unif}}.
  The transformation for the residuals \code{tau} is
  \code{1 - exp(-diff(c(0,tau)))} (cf. Ogata, 1988).
  Another plot inspects the serial correlation between the transformed
  residuals (scatterplot between u_i and u_{i+1}).
}

\usage{
checkResidualProcess(object, plot = 1:2, mfrow = c(1,length(plot)), ...)
}

\arguments{
  \item{object}{
    an object of class \code{"\link{twinSIR}"} or \code{"\link{twinstim}"}.
}
  \item{plot}{
    logical (or integer index) vector indicating if (which) plots of the
    transformed residuals should be produced. The \code{plot} index 1
    corresponds to a \code{\link{ks.plot.unif}} to check for deviations
    of the transformed residuals from the uniform distribution. The
    \code{plot} index 2 corresponds to a scatterplot of \eqn{u_i} vs.
    \eqn{u_{i+1}}. By default (\code{plot = 1:2}), both plots are produced.
  }
  \item{mfrow}{
    see \code{\link{par}}.
  }
  \item{\dots}{
    further arguments passed to \code{\link{ks.plot.unif}}.
  }
}

\value{
  A list (returned invisibly, if \code{plot = TRUE}) with the following
  components:
  \describe{
    \item{tau}{the residual process obtained by
      \code{residuals(object)}.}
    \item{U}{the transformed residuals which should be distributed as
      U(0,1).}
    \item{ks}{the result of the \code{ks.test} for the uniform
      distribution of \code{U}.}
    }
}

\references{
  Ogata, Y. (1988)
  Statistical models for earthquake occurrences and residual analysis
  for point processes.
  \emph{Journal of the American Statistical Association}, 83, 9-27
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{ks.plot.unif}} and the
  \code{\link[=residuals.twinstim]{residuals}}-method for classes
  \code{"twinSIR"} and \code{"twinstim"}.
}

\examples{
data("hagelloch")
fit <- twinSIR(~ household, data = hagelloch)  # a simplistic model
## extract the "residual process", i.e., the fitted cumulative intensities
residuals(fit)
## assess goodness of fit based on these residuals
checkResidualProcess(fit)  # could be better
}

\keyword{dplot}
\keyword{htest}
