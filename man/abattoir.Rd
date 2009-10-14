\name{abattoir}
\alias{abattoir}
\encoding{latin1}

\title{Abattoir Data}
\description{
  A synthetic dataset from the Danish meat inspection -- useful
  for illustrating the beta-binomial CUSUM.
}
\usage{data(abattoir)
}
\details{
  The object of class \code{sts} contains an artificial data set
  inspired by meat inspection data used by Danish Pig Production,
  Denmark. For each week the number of pigs with positive audit reports
  is recorded together with the total number of audits made that week.
}
\seealso{\code{\link{categoricalCUSUM}}}
\examples{
data("abattoir")
plot(abattoir,legend.opts=NULL)
population(abattoir) 
}
\references{
  Höhle, M. (2010), Changepoint detection in categorical time series, 
  Book chapter to appear in T. Kneib and G. Tutz (Eds.), Statistical
  Modelling and Regression Structures, Springer.
}
\keyword{datasets}
