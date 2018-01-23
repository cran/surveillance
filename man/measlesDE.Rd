\name{measlesDE}
\alias{measlesDE}
\docType{data}
\title{Measles in the 16 states of Germany}
\description{
  Weekly number of measles cases in the 16 states
  (Bundeslaender) of Germany for years 2005 to 2007.
}
\usage{data(measlesDE)}
\format{
  An \code{sts} object containing \eqn{156\times 16}{156 x 16}
  observations starting from week 1 in 2005.
  
  The \code{population} slot contains the population fractions
  of each state at 31.12.2006, obtained from the Federal Statistical
  Office of Germany.
}
\source{
  Robert Koch-Institut: SurvStat:
  \url{https://survstat.rki.de/};
  Queried on 14 October 2009. 
}
\seealso{\code{\link{MMRcoverageDE}}}
\examples{
data(measlesDE)
plot(measlesDE)    
}
\references{
  Herzog, S.A., Paul, M. and Held, L. (2011) Heterogeneity in vaccination
    coverage explains the size and occurrence of measles epidemics in
    German surveillance data. Epidemiology and Infection, \bold{139}, 505--515.
}
\keyword{datasets}
