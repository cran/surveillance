\name{fluBYBW}
\alias{fluBYBW}
\encoding{latin1}

\title{Influenza in Southern Germany}
\description{Weekly number of influenza A & B cases in the 140 districts
  of the two Southern German states Bavaria and Baden-Wuerttemberg,
  for years 2001 to 2008.
}
\usage{data(fluBYBW)
}
\format{
  An \code{sts} object containing \eqn{417\times 140}{417 x 140}
  observations starting from week 1 in 2001.
  
  The \code{population} slot contains the population fractions
  of each district at 31.12.2001, obtained from the Federal Statistical
  Office of Germany.
  
  The \code{map} slot contains an object of class \code{SpatialPolygonsDataFrame}.
}
\source{Robert Koch-Institut: SurvStat:
	\url{http://www3.rki.de/SurvStat};
	Queried on 6 March 2009. 
        }
\examples{
    data(fluBYBW)
    # Plot a map of the total number of cases in the year 2001
    plot(fluBYBW[year(fluBYBW) == 2001, ], type = observed ~ 1|unit, 
         labels = FALSE)    
}
\references{
  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136.
}

\keyword{datasets}

