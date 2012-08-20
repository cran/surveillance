\encoding{latin1}
\name{ks.plot.unif}
\alias{ks.plot.unif}
\title{
Plot the ECDF of a uniform sample with Kolmogorov-Smirnov bounds
}

\description{
This plot function takes a univariate sample that should be tested for
a U(0,1) distribution, plots its empirical cumulative distribution
function (\code{\link{ecdf}}), and adds a confidence band by inverting
the corresponding Kolmogorov-Smirnov test (\code{\link{ks.test}}). The
uniform distribution is rejected if the ECDF is not completely inside
the confidence band.
}

\usage{
ks.plot.unif(U, conf.level = 0.95, exact = NULL,
             col.conf = "gray", col.ref = "gray",
             xlab = expression(u[(i)]), ylab = "Cumulative distribution")
}

\arguments{
  \item{U}{
    numeric vector containing the sample.
    Missing values are (silently) ignored.
}
  \item{conf.level}{
    confidence level for the K-S-test (defaults to 0.95), can also be a
    vector of multiple levels.
  }
  \item{exact}{see \code{\link{ks.test}}.}
  \item{col.conf}{
    colour of the confidence lines.
}
  \item{col.ref}{
    colour of the diagonal reference line.
}
  \item{xlab, ylab}{
    axis labels.
  }
}

\value{
\code{NULL} (invisibly).
}

\author{
  Michael H\enc{ö}{oe}hle and Sebastian Meyer.
  The code contains segments originating from the source of the \link{ks.test} function
  \url{http://svn.r-project.org/R/trunk/src/library/stats/R/ks.test.R},
  which is Copyright (C) 1995-2012 The R Core Team available under GPL-2
  (or later) and C functionality from
  \url{http://svn.r-project.org/R/trunk/src/library/stats/src/ks.c},
  which is copyright (C) 1999-2009 the R Core Team and available under
  GPL-2 (or later). Somewhat hidden in their \file{ks.c} file is a statement
  that part of their code is based on code published in George Marsaglia and
  Wai Wan Tsang and Jingbo Wang (2003), "Evaluating Kolmogorov's
  distribution".  Journal of Statistical Software, Volume 8, 2003, Issue 18.
  URL: \url{http://www.jstatsoft.org/v08/i18/}.
}

\seealso{
\code{\link{ks.test}} for the Kolmogorov-Smirnov test, as well as
\code{\link{checkResidualProcess}}, which makes use of this plot
function.
}

\examples{
samp <- runif(99)
ks.plot.unif(samp, conf.level=c(0.95, 0.99), exact=TRUE)
ks.plot.unif(samp, conf.level=c(0.95, 0.99), exact=FALSE)
}

\keyword{hplot}
\keyword{htest}
