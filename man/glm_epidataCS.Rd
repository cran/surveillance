\name{glm_epidataCS}
\alias{glm_epidataCS}
\title{
  Fit an Endemic-Only \code{twinstim} as a Poisson-\code{glm}
}

\description{
  An endemic-only \code{\link{twinstim}} is equivalent to a Poisson
  regression model for the aggregated number of events,
  \eqn{Y_{[t][\bm{s}],k}}, by time-space-type cell. The rate of the
  corresponding Poisson distribution is
  \eqn{e_{[t][\bm{s}]} \cdot \lambda([t],[\bm{s}],k)},
  where \eqn{e_{[t][\bm{s}]} = |[t]| |[\bm{s}]|} is a multiplicative
  offset. Thus, the \code{\link{glm}} function can be used to fit
  an endemic-only \code{twinstim}. However, wrapping in \code{glm} is
  usually slower. 
}

\usage{
glm_epidataCS(formula, data, ...)
}

\arguments{
  \item{formula}{
    an endemic model formula without response, comprising variables of
    \code{data$stgrid} and possibly the variable \code{type} for a
    type-specific model.
  }
  \item{data}{
    an object of class \code{"\link{epidataCS}"}.
  }
  \item{\dots}{
    arguments passed to \code{\link{glm}}. Note that \code{family} and
    \code{offset} are fixed internally.
  }
}

\value{
  a \code{\link{glm}}
}

\author{
  Sebastian Meyer
}

\examples{
data("imdepi")
data("imdepifit")

## Fit an endemic-only twinstim() and an equivalent model wrapped in glm()
fit_twinstim <- update(imdepifit, epidemic = ~0, siaf = NULL,
                       optim.args=list(control=list(trace=0)), verbose=FALSE)
fit_glm <- glm_epidataCS(formula(fit_twinstim)$endemic, imdepi)

## Compare the coefficients
cbind(twinstim=coef(fit_twinstim), glm=coef(fit_glm))
stopifnot(isTRUE(all.equal(coef(fit_glm), coef(fit_twinstim),
                           tolerance = 0.0005, check.attributes = FALSE)))
\dontshow{
if (surveillance.options("allExamples")) {
    ## also check type-specific model:
    stopifnot(isTRUE(all.equal(
        coef(glm_epidataCS(~0+type, imdepi)),
        coef(update(fit_twinstim, endemic=~(1|type))),
    tolerance = 0.0005, check.attributes = FALSE)))
}
}

### also compare to an equivalent endemic-only hhh4() fit

## first need to aggregate imdepi into an "sts" object
load(system.file("shapes", "districtsD.RData", package="surveillance"))
imdepi_sts <- epidataCS2sts(imdepi, freq=12, start=c(2002,1),
    neighbourhood=NULL, tiles=districtsD, popcol.stgrid="popdensity")

## determine the correct offset to get an equivalent model
offset <- 2 * rep(with(subset(imdepi$stgrid, !duplicated(BLOCK)),
                  stop-start), ncol(imdepi_sts)) *
          sum(districtsD$POPULATION) * population(imdepi_sts)

## fit the model using hhh4()
fit_hhh4 <- hhh4(imdepi_sts, control = list(
    end = list(
        f = addSeason2formula(~I(start/365-3.5), period=365, timevar="start"),
        offset = offset
    ), family = "Poisson", subset = 1:nrow(imdepi_sts),
    data = list(start=with(subset(imdepi$stgrid, !duplicated(BLOCK)), start))))

summary(fit_hhh4)
stopifnot(isTRUE(all.equal(coef(fit_hhh4), coef(fit_glm),
                           check.attributes=FALSE)))
}

\keyword{models}
