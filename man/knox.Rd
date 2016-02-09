\name{knox}
\alias{knox}
\alias{plot.knox}
\title{
Knox Test for Space-Time Interaction
}

\description{
Given temporal and spatial distances as well as corresponding critical
thresholds defining what \dQuote{close} means, the function
\code{knox} performs Knox (1963, 1964) test for space-time interaction.
The corresponding p-value can be calculated either by the Poisson
approximation or by a Monte Carlo permutation approach (Mantel, 1967)
with support for parallel computation via \code{\link{plapply}}.
There is a simple \code{plot}-method showing a \code{\link{truehist}} of
the simulated null distribution together with the expected and observed
values. This implementation is due to Meyer et al. (2015).
}

\usage{
knox(dt, ds, eps.t, eps.s, simulate.p.value = TRUE, B = 999, ...)

\method{plot}{knox}(x, ...)
}

\arguments{
  \item{dt,ds}{
    numeric vectors of length \eqn{n*(n-1)/2} containing the pairwise
    temporal and spatial distances of \eqn{n} events.
    Logical vectors indicating temporal/spatial closeness may also be
    supplied, in which case \code{eps.t}/\code{eps.s} is ignored.
  }
  \item{eps.t,eps.s}{
    Critical distances defining closeness in time and space,
    respectively. Distances lower than or equal to the critical distance
    are considered \dQuote{"close"}.
  }
  \item{simulate.p.value}{
    logical indicating if a Monte Carlo permutation test should be
    performed (as per default). Do not forget to set the
    \code{\link{.Random.seed}} via an extra \code{.seed} argument if
    reproducibility is required (see the \dots arguments below).
    If \code{simulate.p.value = FALSE}, the Poisson approximation is
    used (but see the note below).
  }
  \item{B}{
    number of permutations for the Monte Carlo approach.
  }
  \item{\dots}{
    arguments configuring \code{\link{plapply}}:
    \code{.parallel}, \code{.seed}, and \code{.verbose}.
    By default, no parallelization is performed (\code{.parallel = 1}),
    and a progress bar is shown (\code{.verbose = TRUE}).\cr
    For the \code{plot}-method, further arguments passed to
    \code{\link{truehist}}.
  }
  \item{x}{
    an object of class \code{"knox"} as returned by the
    \code{knox} test.
  }
}

\note{
  The Poisson approximation works well if the proportions of close
  pairs in both time and space are small (Kulldorf and Hjalmars,
  1999), otherwise the Monte Carlo permutation approach is
  recommended.
}

\value{
  an object of class \code{"knox"} (inheriting from \code{"htest"}),
  which is a list with the following components:
  \item{method}{a character string indicating the type of test
    performed, and whether the Poisson approximation or Monte Carlo
    simulation was used.}
  \item{data.name}{a character string giving the supplied \code{dt} and
    \code{ds} arguments.}
  \item{statistic}{the number of close pairs.}
  \item{parameter}{if \code{simulate.p.value = TRUE}, the number
    \code{B} of permutations, otherwise the \code{lambda} parameter of
    the Poisson distribution, i.e., the same as \code{null.value}.}
  \item{p.value}{the p-value for the test. In case
    \code{simulate.p.value = TRUE}, the p-value from the Poisson
    approximation is still attached as an attribute \code{"Poisson"}.}
  \item{alternative}{the character string \code{"greater"} (this is a
    one-sided test).}
  \item{null.value}{the expected number of close pairs in the absence of
    space-time interaction.}
  \item{table}{the contingency table of \code{dt <= eps.t} and
    \code{ds <= eps.s}.}

  The \code{plot}-method invisibly returns \code{NULL}.
}

\author{
  Sebastian Meyer
}

\seealso{
  the space-time K-function test \code{\link{stKtest}}
  to combine information from different scales 
  \code{eps.t} and \code{eps.s} while also handling edge effects,
  as well as function \code{\link{epitest}} for testing
  \code{"\link{twinstim}"} models.
}

\references{
  Knox, G. (1963):
  Detection of low intensity epidemicity: application to cleft lip and palate.
  \emph{British Journal of Preventive & Social Medicine}, \bold{17}, 121-127.
  
  Knox, E. G. (1964):
  The detection of space-time interactions.
  \emph{Journal of the Royal Statistical Society. Series C (Applied
    Statistics)}, \bold{13}, 25-30.

  Kulldorff, M. and Hjalmars, U. (1999):
  The Knox method and other tests for space-time interaction.
  \emph{Biometrics}, \bold{55}, 544-552.
  
  Mantel, N. (1967):
  The detection of disease clustering and a generalized regression approach.
  \emph{Cancer Research}, \bold{27}, 209-220.

  Meyer, S., Warnke, I., R\enc{ö}{oe}ssler, W. and Held, L. (2015):
  Model-based testing for space-time interaction using point processes:
  An application to psychiatric hospital admissions in an urban area.
  \emph{Spatial and Spatio-temporal Epidemiology}. In revision.\cr
  Preprint available at \url{http://arxiv.org/abs/1512.09052}
}

\examples{
data("imdepi")
imdepiB <- subset(imdepi, type == "B")

## Obtain the p-value via a Monte Carlo permutation test,
## where the permutations can be computed in parallel
## (using forking on Unix-alikes and a cluster on Windows, see ?plapply)
knoxtest <- knox(
    dt = dist(imdepiB$events$time), eps.t = 30,
    ds = dist(coordinates(imdepiB$events)), eps.s = 50,
    simulate.p.value = TRUE, B = 199,
    .parallel = 2, .seed = 1, .verbose = FALSE
)
knoxtest
plot(knoxtest)
}

\keyword{htest}
