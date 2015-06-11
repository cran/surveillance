\name{twinstim_epitest}
\alias{epitest}
\alias{coef.epitest}
\alias{plot.epitest}
\title{Permutation Test for Space-Time Interaction in \code{"twinstim"}}

\description{
  The function \code{epitest} takes a simple epidemic \code{"twinstim"}
  model (one with \code{epidemic = ~1}) and tests if the spatio-temporal
  interaction invoked by the \code{epidemic} model component is
  statistically significant.
}

\usage{
epitest(model, data, tiles, method = "time", B = 199,
        eps.s = NULL, eps.t = NULL, fixed = NULL,
        verbose = TRUE, compress = FALSE, ...)

\method{coef}{epitest}(object, which = c("m1", "m0"), ...)
\method{plot}{epitest}(x, teststat = c("simpleR0", "D"), ...)
}

\arguments{
  \item{model}{
    a simple epidemic \code{"\link{twinstim}"} without covariates, i.e.,
    \code{epidemic = ~1}. This is because covariate effects in the
    epidemic component are not well identified when there is no
    space-time interaction such as in the permuted data. Estimating a
    rich epidemic \code{model} under the null hypothesis of no
    space-time interaction will most likely result in singular convergence.
  }
  \item{data}{
    an object of class \code{"\link{epidataCS}"}, the \code{data} to
    which the \code{model} was fitted.
  }
  \item{tiles}{
    (only used by \code{method = "simulate"})
    a \code{"\linkS4class{SpatialPolygons}"} representation of the
    \code{tile}s in \code{data$stgrid}.
  }
  \item{method}{
    one of the following character strings specifying the test method:
    \describe{
      \item{\code{"LRT"}:}{
        a simple likelihood ratio test of the epidemic
        \code{model} against the corresponding endemic-only model,
      }
      \item{\code{"time"}/\code{"space"}:}{
        a Monte Carlo permutation test where the null distribution is
        obtained by relabeling time points or locations, respectively
        (using \code{\link{permute.epidataCS}}).
      }
      \item{\code{"simulate"}:}{
        obtain the null distribution of the test statistic by
        simulations from the endemic-only model
        (using \code{\link{simulate.twinstim}}).
      }
    }
  }
  \item{B}{
    the number of permutations for the Monte Carlo approach.
    The default number is rather low; if computationally feasible,
    \code{B = 999} is more appropriate. Note that this determines the
    \dQuote{resolution} of the p-value: the smallest attainable p-value
    is \code{1/(B+1)}.
  }
  \item{eps.s,eps.t}{arguments for \code{\link{simpleR0}}.}
  \item{fixed}{
    optional character vector naming parameters to fix at their original
    value when re-fitting the \code{model} on permuted data.
    The special value \code{fixed = TRUE} means to fix all epidemic
    parameters but the intercept.
  }
  \item{verbose}{
    the amount of tracing in the range \code{0:3}.
    Set to 0 (or \code{FALSE}) for no output,
    1 (or \code{TRUE}, the default) for a progress bar,
    2 for the test statistics resulting from each permutation,
    and to 3 for additional tracing of the log-likelihood
    maximization in each permutation (not useful if parallelized).
    Tracing does not work if permutations are parallelized using clusters.
    See \code{\link{plapply}} for other choices.
  }
  \item{compress}{
    logical indicating if the \code{nobs}-dependent elements \code{"fitted"},
    \code{"fittedComponents"}, and \code{"R0"} should be dropped from
    the permutation-based model fits. Not keeping these elements saves a
    lot of memory especially with a large number of events.
    Note, however, that the returned \code{permfits} then no longer are
    fully valid \code{"twinstim"} objects (but most methods will still work).
  }
  \item{\dots}{further arguments for \code{\link{plapply}} to configure
    parallel operation, i.e., \code{.parallel} as well as
    \code{.seed} to make the results reproducible.\cr
    For the \code{plot}-method, further arguments passed to
    \code{\link{truehist}}.\cr
    Ignored by the \code{coef}-method.
  }
  \item{object,x}{
    an object of class \code{"epitest"} as returned by \code{epitest}.
  }
  \item{which}{
    a character string indicating either the full (\code{"m1"}, default)
    or the endemic-only (\code{"m0"}) model.
  }
  \item{teststat}{
    a character string determining the test statistic to plot, either
    \code{"\link{simpleR0}"} or \code{"D"} (twice the log-likelihood
    difference of the models).
  }
}

\value{
  a list (inheriting from \code{"htest"}) with the following components:
  
  \item{method}{a character string indicating the type of test performed.}
  \item{data.name}{a character string giving the supplied \code{data} and
    \code{model} arguments.}
  \item{statistic}{the observed test statistic.}
  \item{parameter}{the (effective) number of permutations used to
    calculate the p-value (only those with convergent fits are used).}
  \item{p.value}{the p-value for the test. For the \code{method}s
    involving resampling under the null (\code{method != "LRT"}),
    it is based on the subset of convergent fits only and the p-value
    from the simple LRT is attached as an attribute \code{"LRT"}.}

  In addition, if \code{method != "LRT"}, the result will have the
  following elements:
  
  \item{permfits}{the list of model fits (endemic-only and epidemic)
    from the \code{B} permutations.}
  \item{permstats}{a data frame with \code{B} rows and the columns
    \code{"l0"} (log-likelihood of the endemic-only model \code{m0}),
    \code{"l1"} (log-likelihood of the epidemic model \code{m1}),
    \code{"D"} (twice their difference),
    \code{"simpleR0"} (the results of \code{\link{simpleR0}(m1, eps.s, eps.t)}), 
    and \code{"converged"} (a boolean indicator if both models converged).}

  The \code{plot}-method invisibly returns \code{NULL}.
  The \code{coef}-method returns the \code{B} x \code{length(coef(model))}
  matrix of parameter estimates.
}

\details{
  The test statistic is the reproduction number \code{\link{simpleR0}}.
  A likelihood ratio test of the supplied epidemic model against
  the corresponding endemic-only model is also available.
  The null distribution of the test statistic is obtained by a Monte Carlo permutation
  approach (via \code{\link{permute.epidataCS}}) similar to the
  \code{\link{knox}} test.
  
  The \code{plot}-method shows a \code{\link{truehist}} of
  the simulated null distribution together with the observed value.
  The \code{coef}-method extracts the parameter estimates from the \code{B}
  \code{permfits} (by default for the full model \code{which = "m1"}).
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{permute.epidataCS}}, \code{\link{knox}}
}

\examples{
data("imdepi")
data("imdepifit")

## test for space-time interaction of the B-cases
## assuming spatial interaction to be constant within 50 km
imdepiB50 <- update(subset(imdepi, type == "B"), eps.s = 50)
imdfitB50 <- update(imdepifit, data = imdepiB50,
                    epidemic = ~1, epilink = "identity",
                    siaf = NULL, control.siaf = NULL,
                    start = c("e.(Intercept)" = 1e-6))

## simple likelihood ratio test
epitest(imdfitB50, imdepiB50, method = "LRT")

## permutation test (only few permutations for speed, in parallel)
et <- epitest(imdfitB50, imdepiB50, B = 4 + 25*surveillance.options("allExamples"),
              verbose = 2 * (.Platform$OS.type == "unix"),
              .seed = 1, .parallel = 2)
et
plot(et)
summary(coef(et))
}

\keyword{htest}
