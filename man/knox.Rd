\encoding{latin1}
\name{knox}
\alias{knox}
\alias{plot.knox}
\alias{toLatex.knox}
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
values.
This implementation of the Knox test is due to Meyer et al. (2016).
}

\usage{
knox(dt, ds, eps.t, eps.s, simulate.p.value = TRUE, B = 999, ...)

\method{plot}{knox}(x, ...)
}

\arguments{
  \item{dt,ds}{
    numeric vectors containing temporal and spatial distances, respectively.
    Logical vectors indicating temporal/spatial closeness may also be
    supplied, in which case \code{eps.t}/\code{eps.s} is ignored.
    To test for space-time interaction in a single point pattern of \eqn{n}
    events, these vectors should be of length \eqn{n*(n-1)/2} and
    contain the pairwise event distances (e.g., the lower triangle of
    the distance matrix, such as in \code{"\link{dist}"} objects).
    Note that there is no special handling of matrix input, i.e.,
    if \code{dt} or \code{ds} are matrices, all elements are used
    (but a warning is given if a symmetric matrix is detected).
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
  pairs in both time and space are small (Kulldorff and Hjalmars,
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

  A \code{toLatex}-method exists, which generates LaTeX code for the
  contingency table associated with the Knox test.
}

\author{
  Sebastian Meyer
}

\seealso{
  The function \code{mantel.randtest} in package \pkg{ade4}
  implements Mantel's (1967) space-time interaction test, i.e., using
  the Pearson correlation between the spatial and temporal distances of
  all event pairs as the test statistic, and assessing statistical
  significance using a Monte Carlo permutation approach as with
  \code{simulate.p.value} here in the \code{knox} function.
  To combine information from different scales \code{eps.t} and
  \code{eps.s} while also handling edge effects, the space-time
  K-function test available via \code{\link{stKtest}} can be used.
  Function \code{\link{epitest}} tests epidemicity in a
  \code{"\link{twinstim}"} point process model.
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

  Meyer, S., Warnke, I., R\enc{ö}{oe}ssler, W. and Held, L. (2016):
  Model-based testing for space-time interaction using point processes:
  An application to psychiatric hospital admissions in an urban area.
  \emph{Spatial and Spatio-temporal Epidemiology}, \bold{17}, 15-25.
  \doi{10.1016/j.sste.2016.03.002}.
  Eprint: \url{http://arxiv.org/abs/1512.09052}.
}

\examples{
data("imdepi")
imdepiB <- subset(imdepi, type == "B")

## Perfom the Knox test using the Poisson approximation
knoxtest <- knox(
    dt = dist(imdepiB$events$time), eps.t = 30,
    ds = dist(coordinates(imdepiB$events)), eps.s = 50,
    simulate.p.value = FALSE
)
knoxtest
## The Poisson approximation works well for these data since
## the proportion of close pairs is rather small (204/56280).

## contingency table in LaTeX
toLatex(knoxtest)


if (surveillance.options("allExamples")) {
    ## Obtain the p-value via a Monte Carlo permutation test,
    ## where the permutations can be computed in parallel
    ## (using forking on Unix-alikes and a cluster on Windows, see ?plapply)
    knoxtestMC <- knox(
        dt = dist(imdepiB$events$time), eps.t = 30,
        ds = dist(coordinates(imdepiB$events)), eps.s = 50,
        simulate.p.value = TRUE, B = 999,
        .parallel = 2, .seed = 1, .verbose = FALSE
    )
    knoxtestMC
    plot(knoxtestMC)
}
}

\keyword{htest}
