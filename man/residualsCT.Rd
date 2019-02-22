\name{residualsCT}
\alias{residuals.twinSIR}
\alias{residuals.twinstim}
\alias{residuals.simEpidataCS}

\title{
Extract Cox-Snell-like Residuals of a Fitted Point Process
}

\description{
Extract the \dQuote{residual process} (cf. Ogata, 1988) of a fitted
point process model specified through the conditional intensity
function, for instance a model of class \code{"\link{twinSIR}"} or
\code{"\link{twinstim}"} (and also \code{"\link{simEpidataCS}"}).
The residuals are defined as the fitted cumulative intensities at the
event times, and are generalized residuals similar to those discussed in
Cox and Snell (1968).
}

\usage{
\method{residuals}{twinSIR}(object, ...)
\method{residuals}{twinstim}(object, ...)
\method{residuals}{simEpidataCS}(object, ...)
}

\arguments{
  \item{object}{
    an object of one of the aforementioned model classes.
  }
  \item{\dots}{unused (argument of the generic).}
}

\details{
  For objects of class \code{twinstim}, the residuals may already be
  stored in the object as component \code{object$tau} if the model was
  fitted with \code{cumCIF = TRUE} (and they always are for
  \code{"simEpidataCS"}). In this case, the \code{residuals}
  method just extracts these values.  Otherwise, the residuals have to
  be calculated, which is only possible with access to the model
  environment, i.e. \code{object} must have been fitted with
  \code{model = TRUE}. The calculated residuals are then also appended
  to \code{object} for future use. However, if \code{cumCIF} and
  \code{model} were both set to true in the \code{object} fit, then it
  is not possible to calculate the residuals and the method returns an
  error.
}

\value{
Numeric vector of length the number of events of the corresponding point
process fitted by \code{object}. This is the observed residual process.
}

\references{
  Ogata, Y. (1988)
  Statistical models for earthquake occurrences and residual analysis
  for point processes.
  \emph{Journal of the American Statistical Association}, 83, 9-27

  Cox, D. R. & Snell, E. J. (1968)
  A general definition of residuals.
  \emph{Journal of the Royal Statistical Society. Series B (Methodological)}, 30, 248-275
}

\seealso{
  \code{\link{checkResidualProcess}} to graphically check the
  goodness-of-fit of the underlying model.
}

\author{
  Sebastian Meyer
}

\keyword{methods}
