\name{hhh4_predict}
\alias{predict.hhh4}
\title{Predictions from a \code{hhh4} Model}
\description{
  Get fitted (component) means from a \code{\link{hhh4}} model.
}

\usage{
\method{predict}{hhh4}(object, newSubset=object$control$subset,
        type="response", \dots)
}

\arguments{
\item{object}{fitted \code{\link{hhh4}} model (class \code{"hhh4"}).}
\item{newSubset}{subset of time points for which to return the
  predictions. Defaults to the subset used for fitting the model, and
  must be a subset of \code{1:nrow(object$stsObj)}.}
\item{type}{the type of prediction required. The default
  (\code{"response"} or, equivalently, \code{"mean"}) is on the 
  scale of the response variable (mean = endemic plus epidemic components).
  The alternatives are: \code{"endemic"}, \code{"epidemic"},
  \code{"epi.own"} (i.e. the autoregressive part), and
  \code{"epi.neighbours"} (i.e. the spatio-temporal part).}
\item{\dots}{unused (argument of the generic).}
}

\value{
  matrix of fitted means for each time point (of \code{newSubset}) and region.
}

\note{
  Predictions for \dQuote{newdata}, i.e., with modified covariates or
  fixed weights, can be computed manually by adjusting the control list
  (in a copy of the original fit), dropping the old \code{terms}, and using
  the internal function \code{\link{meanHHH}} directly, see the Example.
}

\author{Michaela Paul and Sebastian Meyer}

\keyword{methods}
\keyword{models}

\examples{
## simulate simple seasonal noise with reduced baseline for t >= 60
t <- 0:100
y <- rpois(length(t), exp(3 + sin(2*pi*t/52) - 2*(t >= 60)))
obj <- sts(y)
plot(obj)

## fit true model
fit <- hhh4(obj, list(end = list(f = addSeason2formula(~lock)),
                      data = list(lock = as.integer(t >= 60)),
                      family = "Poisson"))
coef(fit, amplitudeShift = TRUE, se = TRUE)

## compute predictions for a subset of the time points
stopifnot(identical(predict(fit), fitted(fit)))
plot(obj)
lines(40:80, predict(fit, newSubset = 40:80), lwd = 2)


## advanced: compute predictions for "newdata" (here, a modified covariate)
mod <- fit
mod$terms <- NULL  # to be sure
mod$control$data$lock[t >= 60] <- 0.5
pred <- meanHHH(mod$coefficients, terms(mod))$mean
plot(fit, xaxis = NA)
lines(mod$control$subset, pred, lty = 2)
}
