\name{fluBYBW}
\alias{fluBYBW}
\docType{data}

\title{Influenza in Southern Germany}

\description{
  Weekly number of influenza A & B cases in the 140 districts
  of the two Southern German states Bavaria and Baden-Wuerttemberg,
  for the years 2001 to 2008. These surveillance data have been
  analyzed originally by Paul and Held (2011) and more recently by
  Meyer and Held (2014).
}

\usage{data(fluBYBW)}

\format{
  An \code{sts} object containing \eqn{416\times 140}{416 x 140}
  observations starting from week 1 in 2001.
  
  The \code{population} slot contains the population fractions
  of each district at 31.12.2001, obtained from the Federal Statistical
  Office of Germany.
  
  The \code{map} slot contains an object of class
  \code{"\linkS4class{SpatialPolygonsDataFrame}"}.
}

\source{
  Robert Koch-Institut: SurvStat: \url{https://survstat.rki.de/};
  Queried on 6 March 2009.
}

\note{
  Prior to \pkg{surveillance} version 1.6-0, \code{data(fluBYBW)}
  contained a redundant last row (417) filled with zeroes only.
}

\examples{
data("fluBYBW")
# Count time series plot
plot(fluBYBW, type = observed ~ time)
# Map of disease incidence (per 100000 inhabitants) for the year 2001
plot(fluBYBW, type = observed ~ unit, tps = 1:52, total.args = list(),
     population = fluBYBW@map$X31_12_01 / 100000)
# the overall rate for 2001 shown in the bottom right corner is
sum(observed(fluBYBW[1:52,])) / sum(fluBYBW@map$X31_12_01) * 100000

\dontrun{
# Generating an animation takes a while.
# Here we take the first 20 weeks of 2001 (runtime: ~3 minutes).
# The full animation is available in Supplement A of Meyer and Held (2014)
if (require("animation")) {
    oldwd <- setwd(tempdir())  # to not clutter up the current working dir
    saveHTML(animate(fluBYBW, tps = 1:20),
             title="Evolution of influenza in Bayern and Baden-Wuerttemberg",
             ani.width=500, ani.height=600)
    setwd(oldwd)
}
}
}

\references{
  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118-1136.

  Meyer, S. and Held, L. (2014):
    Power-law models for infectious disease spread.
    \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.
    \doi{10.1214/14-AOAS743}
}

\keyword{datasets}
