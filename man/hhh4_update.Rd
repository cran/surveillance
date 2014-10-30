\name{hhh4_update}
\alias{update.hhh4}
\title{
\code{update} a fitted \code{"hhh4"} model
}

\description{
  Re-fit a \code{"\link{hhh4}"} model with a modified \code{control} list.
}

\usage{
\method{update}{hhh4}(object, ..., S = NULL, subset.upper = NULL,
       use.estimates = FALSE, evaluate = TRUE)
}

\arguments{
  \item{object}{
    a fitted \code{"hhh4"} model.
  }
  \item{\dots}{
    components modifying the original control list for
    \code{\link{hhh4}}. Modifications are performed by
    \code{\link{modifyList}(object$control, list(...))}.
  }
  \item{S}{
    a named list of numeric vectors serving as argument for
    \code{\link{addSeason2formula}}, or \code{NULL} (meaning no
    modification of seasonal terms). This argument provides a convenient
    way of changing the number of harmonics in the \code{f}ormulae of
    the model components \code{"ar"}, \code{"ne"} and \code{"end"} (to be
    used as names of the list). Non-specified components are not touched.
    Updating the \code{i}'th component's \code{f}ormula works by first
    dropping all sine and cosine terms and then applying
    \code{addSeason2formula} with arguments \code{S = S[[i]]} and
    \code{period = object$stsObj@freq}. Note that this step of updating
    seasonality is processed after modification of the \code{control}
    list by the \code{\dots} arguments.
  }
  \item{subset.upper}{
    if a scalar value, restrict the new fit to those epochs of
    \code{object$control$subset} which are \code{<= subset.upper}.
    This argument is used by \code{\link{oneStepAhead}}.
  }
  \item{use.estimates}{
    logical specifying if the \code{coef(object)} should be used as
    starting values for the new fit. This currently only works if the
    set of parameters is unchanged by the updated model formulation and
    thus defaults to \code{FALSE}.
  }
  \item{evaluate}{
    logical indicating if the updated model should be fitted directly
    (defaults to \code{TRUE}). Otherwise, the updated \code{control} list
    is returned.
  }
}

\value{
  If \code{evaluate = TRUE} the re-fitted object, otherwise the updated
  \code{control} list for \code{\link{hhh4}}.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{hhh4}}
}

\examples{
data("salmonella.agona")
## convert to sts class
salmonella <- disProg2sts(salmonella.agona)

## fit a basic model
fit0 <- hhh4(salmonella,
            list(ar = list(f = ~1), end = list(f = addSeason2formula(~t))))

## update: Poisson -> NegBin1, component seasonality
fit1 <- update(fit0, family = "NegBin1", S = list(end=2, ar=2))

## compare fits
AIC(fit0, fit1)
opar <- par(mfrow=c(2,2))
plot(fit0, type="fitted", names="fit0", par.settings=NULL)
plot(fit1, type="fitted", names="fit1", par.settings=NULL)
plot(fit0, fit1, type="season", components=c("end", "ar"), par.settings=NULL)
par(opar)
}

\keyword{models}
\keyword{methods}
