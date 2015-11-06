\name{hhh4_calibration}
\alias{calibrationTest.hhh4}
\title{Test Calibration of a \code{hhh4} Model}

\description{
  The \code{\link{calibrationTest}}-method for \code{"\link{hhh4}"}
  models tests the calibration of the model fit, i.e., if the observed
  counts are from the fitted Poisson/NegBin distribution.
}

\usage{
\method{calibrationTest}{hhh4}(x,
                subset = x$control$subset, units = seq_len(x$nUnit), ...)
}

\arguments{
  \item{x}{
    an object of class \code{"hhh4"}.
  }
  \item{subset}{
    subset of time points for which to test calibration.
    Defaults to the subset used for fitting the model.
  }
  \item{units}{
    integer or character vector selecting the \code{units} for which to
    test calibration. It indexes \code{colnames(x$stsObj)}.
    By default all units are taken into account.
  }
  \item{\dots}{
    arguments passed to \code{\link{calibrationTest.default}},
    e.g., \code{which} to select a scoring rule.
  }
}

\value{
  an object of class \code{"htest"}, see \code{\link{calibrationTest.default}}.
}

\author{
  Sebastian Meyer
}

\examples{
data("measlesWeserEms")

## fit a hhh4 model
measlesModel <- list(
    end = list(f = addSeason2formula(~0 + ri(type="iid"))),
    ar = list(f = ~1),
    ne = list(f = ~1, weights = W_powerlaw(maxlag = 5)),
    family = "NegBin1")
measlesFit <- hhh4(measlesWeserEms, control = measlesModel)

## plot fitted mean components
plot(measlesFit, units = NULL)

## overall calibration test
calibrationTest(measlesFit) # default is DSS (not suitable for low counts)
calibrationTest(measlesFit, which = "rps")
calibrationTest(measlesFit, which = "logs")
                
## test calibration for the second year in a specific district
calibrationTest(measlesFit, subset = 53:104, units = "03452", which = "rps")
}

\keyword{htest}
\keyword{methods}
