\encoding{latin1}
\name{twinstim_epitest}
\alias{epitest}
\alias{coef.epitest}
\alias{plot.epitest}
\title{Permutation Test for Space-Time Interaction in \code{"twinstim"}}

\description{
  The function \code{epitest} takes a \code{"twinstim"} model
  and tests if the spatio-temporal interaction invoked by the epidemic
  model component is statistically significant.
  The test only works for simple epidemic models, where \code{epidemic = ~1}
  (no additional parameters for event-specific infectivity),
  and requires the non-canonical \code{epilink="identity"} (see
  \code{\link{twinstim}}).
  A permutation test is performed by default, which is only valid if the
  endemic intensity is space-time separable.
  The approach is described in detail in Meyer et al. (2016),
  where it is also compared to alternative global tests for clustering
  such as the \code{\link{knox}} test.
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
    a simple epidemic \code{"\link{twinstim}"} with \code{epidemic = ~1},
    fitted using the non-canonical \code{epilink="identity"}.
    Note that the permutation test is only valid for models with
    a space-time separable endemic intensity, where covariates vary
    either in space or time but not both.
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
        (using \code{\link{simEndemicEvents}}).
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
  This space-time interaction test is limited to models with
  \code{epidemic = ~1}, since covariate effects are not identifiable
  under the null hypothesis of no space-time interaction.
  Estimating a rich epidemic \code{model} based on permuted data
  will most likely result in singular convergence.
  A similar issue might arise when the model employs parametric
  interaction functions, in which case \code{fixed=TRUE} can be used.
  For further details see Meyer et al. (2016).
  
  The test statistic is the reproduction number \code{\link{simpleR0}}.
  A likelihood ratio test of the supplied epidemic model against
  the corresponding endemic-only model is also available.
  By default, the null distribution of the test statistic under no
  space-time interaction is obtained by a Monte Carlo permutation
  approach (via \code{\link{permute.epidataCS}}) and therefore relies on
  a space-time separable endemic model component.
  
  The \code{plot}-method shows a \code{\link{truehist}} of
  the simulated null distribution together with the observed value.
  The \code{coef}-method extracts the parameter estimates from the \code{B}
  \code{permfits} (by default for the full model \code{which = "m1"}).
}

\references{
  Meyer, S., Warnke, I., R\enc{ö}{oe}ssler, W. and Held, L. (2016):
  Model-based testing for space-time interaction using point processes:
  An application to psychiatric hospital admissions in an urban area.
  \emph{Spatial and Spatio-temporal Epidemiology}, \bold{17}, 15-25.
  \doi{10.1016/j.sste.2016.03.002}.
  Eprint: \url{http://arxiv.org/abs/1512.09052}.
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{permute.epidataCS}}, \code{\link{knox}}
}

\examples{
data("imdepi", "imdepifit")

## test for space-time interaction of the B-cases
## assuming spatial interaction to be constant within 50 km
imdepiB50 <- update(subset(imdepi, type == "B"), eps.s = 50)
imdfitB50 <- update(imdepifit, data = imdepiB50,
                    epidemic = ~1, epilink = "identity", siaf = NULL,
                    start = c("e.(Intercept)" = 0))

## simple likelihood ratio test
epitest(imdfitB50, imdepiB50, method = "LRT")

## permutation test (only a few permutations for speed)
et <- epitest(imdfitB50, imdepiB50, B = 3 + 26*surveillance.options("allExamples"),
              verbose = 2 * (.Platform$OS.type == "unix"),
              .seed = 1, .parallel = 1 + surveillance.options("allExamples"))
et
plot(et)
## evidence against the null hypothesis of no space-time interaction

summary(coef(et, which = "m1"))
}

\keyword{htest}
