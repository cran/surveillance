\name{permutationTest}
\alias{permutationTest}
\title{Monte Carlo Permutation Test for Paired Individual Scores}

\description{
  The difference between mean \code{\link{scores}} from model 1 and
  mean \code{\link{scores}} from model 2 is used as the test statistic.
  Under the null hypothesis of no difference, the actually observed
  difference between mean scores should not be notably different from
  the distribution of the test statistic under permutation.
  As the computation of all possible permutations is only feasible for
  small datasets, a random sample of permutations is used to obtain the
  null distribution. The resulting p-value thus depends on the
  \code{\link{.Random.seed}}.
}

\usage{
permutationTest(score1, score2, nPermutation = 9999,
                plot = FALSE, verbose = FALSE)
}

\arguments{
  \item{score1, score2}{
    numeric vectors of scores from models 1 and 2, respectively.
  }
  \item{nPermutation}{
    number of Monte Carlo replicates.
  }
  \item{plot}{
    logical indicating if a \code{\link[MASS]{truehist}} of the \code{nPermutation}
    permutation test statistics should be plotted with a vertical line
    marking the observed difference of the means.
    To customize the histogram, \code{plot} can also be a list of
    arguments for \code{truehist} replacing internal defaults.
  }
  \item{verbose}{
    logical indicating if the results should be printed in one line.
  }
}

\details{
  For each permutation, we first randomly assign the membership of the n
  individual scores to either model 1 or 2 with probability 0.5. We then
  compute the respective difference in mean for model 1 and 2 in this
  permuted set of scores. The Monte Carlo p-value is then given by
  (1 + #\{permuted differences larger than observed difference (in
  absolute value)\}) / (1 + \code{nPermutation}).
}

\value{
  a list of the following elements:
  \item{diffObs}{observed difference in mean scores, i.e.,
    \code{mean(score1) - mean(score2)}}
  \item{pVal.permut}{p-value of the permutation test}
  \item{pVal.t}{p-value of the corresponding
    \code{\link{t.test}(score1, score2, paired=TRUE)}}
}

\author{
  Michaela Paul with contributions by Sebastian Meyer
}

\references{
  Paul, M. and Held, L. (2011):
  Predictive assessment of a non-linear random effects model for
  multivariate time series of infectious disease counts.
  \emph{Statistics in Medicine}, \bold{30} (10), 1118-1136.
  \doi{10.1002/sim.4177}
}

\seealso{
  Package \CRANpkg{coin} for a comprehensive permutation test framework.
}

\examples{
permutationTest(rnorm(50, 1.5), rnorm(50, 1), plot = TRUE)
}

\keyword{htest}
