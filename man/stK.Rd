\encoding{latin1}
\name{stK}
\alias{stKtest}
\alias{plot.stKtest}
\title{
Diggle et al (1995) K-function test for space-time clustering
}

\description{
The function \code{stKtest} wraps functions in package \pkg{splancs} to
perform the K-function based Monte Carlo permutation test for space-time
clustering (Diggle et al, 1995) for \code{"epidataCS"}.
The implementation is due to Meyer et al. (2016).
}

\usage{
stKtest(object, eps.s = NULL, eps.t = NULL, B = 199,
        cores = 1, seed = NULL, poly = object$W)

\method{plot}{stKtest}(x, which = c("D", "R", "MC"),
     args.D = list(), args.D0 = args.D, args.R = list(), args.MC = list(),
     mfrow = sort(n2mfrow(length(which))), ...)
}

\arguments{
  \item{object}{an object of class \code{"epidataCS"}.}
  \item{eps.s, eps.t}{
    numeric vectors defining the spatial and temporal
    grids of critical distances over which to evaluate the test.
    The default (\code{NULL}) uses equidistant values from 0 to the
    smallest \code{eps.s}/\code{eps.t} value in \code{object$events},
    but not larger than half the observed spatial/temporal domain.
  }
  \item{B}{the number of permutations.}
  \item{cores}{
    the number of parallel processes over which to distribute the
    requested number of permutations.
  }
  \item{seed}{
    argument for \code{\link{set.seed}} to initialize the random number
    generator such that results become reproducible
    (also if \code{cores > 1}, see \code{\link{plapply}}).
  }
  \item{poly}{
    the polygonal observation region of the events (as an object handled
    by \code{\link{xylist}}). The default \code{object$W} might not work
    since package \pkg{splancs} does not support multi-polygons. In this
    case, the \code{poly} argument can be used to specify a substitute.
  }
  \item{x}{an \code{"stKtest"}.}
  \item{which}{
    a character vector indicating which diagnostic plots to produce.
    The full set is \code{c("D", "D0", "R", "MC")}.
    The special value \code{which = "stdiagn"} means to call the
    associated \pkg{splancs} function \code{\link[splancs]{stdiagn}}.
  }
  \item{args.D,args.D0,args.R,args.MC}{
    argument lists for the plot functions \code{\link{persp}} (for
    \code{"D"} and \code{"D0"}), \code{\link{plot.default}}
    (\code{"R"}), and \code{\link[MASS]{truehist}} (\code{"MC"}),
    respectively, to modify the default settings.
    Ignored if \code{which = "stdiagn"}.
  }
  \item{mfrow}{
    \code{\link{par}}-setting to layout the plots.
    Ignored for \code{which = "stdiagn"} and if set to \code{NULL}.
  }
  \item{\dots}{ignored (argument of the generic).}
}

\value{
  an object of class \code{"stKtest"} (inheriting from \code{"htest"}),
  which is a list with the following components:
  \item{method}{a character string indicating the type of test performed.}
  \item{data.name}{a character string naming the supplied \code{object}.}
  \item{statistic}{the sum \eqn{U} of the standardized residuals \eqn{R(s,t)}.}
  \item{parameter}{the number \code{B} of permutations.}
  \item{p.value}{the p-value for the test.}
  \item{pts}{the coordinate matrix of the event locations (for
    \code{\link[splancs]{stdiagn}}.}
  \item{stK}{the estimated K-function as returned by
    \code{\link[splancs]{stkhat}}.}
  \item{seD}{the standard error of the estimated \eqn{D(s,t)} as
    returned by \code{\link[splancs]{stsecal}}.}
  \item{mctest}{the observed and permutation values of the test
    statistic as returned by \code{\link[splancs]{stmctest}}.}

  The \code{plot}-method invisibly returns \code{NULL}.
}

\references{
  Diggle, P. J.; Chetwynd, A. G.; H\enc{ä}{ae}ggkvist, R. and Morris, S. E. (1995):
  Second-order analysis of space-time clustering
  \emph{Statistical Methods in Medical Research}, \bold{4}, 124-136.

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
  the simple \code{\link{knox}} test and function \code{\link{epitest}}
  for testing \code{"\link{twinstim}"} models.
}

\examples{
if (requireNamespace("splancs")) {
    data("imdepi")
    imdepiB <- subset(imdepi, type == "B")
    mainpoly <- coordinates(imdepiB$W@polygons[[1]]@Polygons[[5]])
    if (surveillance.options("allExamples")) {
        SGRID <- c(0, 10, 25, 50, 75, 100, 150, 200)
        TGRID <- c(0, 7, 14, 21, 28)
        B <- 99
        CORES <- 2
    } else { # dummy settings for fast CRAN checks
        SGRID <- c(0, 50)
        TGRID <- c(0, 30)
        B <- 9
        CORES <- 1
    }
    imdBstKtest <- stKtest(imdepiB, eps.s = SGRID, eps.t = TGRID, B = B,
                           cores = CORES, seed = 1, poly = list(mainpoly))
    print(imdBstKtest)
    plot(imdBstKtest)
}
}

\keyword{htest}
