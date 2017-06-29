\name{algo.farrington.threshold}
\alias{algo.farrington.threshold}
\title{Compute prediction interval for a new observation}
\description{
Depending on the current transformation \eqn{h(y)= \{y, \sqrt{y}, y^{2/3}\}},

\deqn{V(h(y_0)-h(\mu_0))=V(h(y_0))+V(h(\mu_0))}

is used to compute a prediction interval. The prediction variance
consists of a component due to the variance of having a single
observation and a prediction variance.
}
\usage{
algo.farrington.threshold(pred,phi,alpha=0.01,skewness.transform="none",y)
}
\arguments{
\item{pred}{A GLM prediction object}
\item{phi}{Current overdispersion parameter (superflous?)}
\item{alpha}{Quantile level in Gaussian based CI, i.e. an \eqn{(1-\alpha)\cdot 100\%}
    confidence interval is computed. }
\item{skewness.transform}{Skewness correction, i.e. one of
    \code{"none"}, \code{"1/2"}, or \code{"2/3"}.}
\item{y}{Observed number}
}
\value{
  Vector of length four with lower and upper bounds of an
  \eqn{(1-\alpha)\cdot 100\%} confidence interval (first two
  arguments) and corresponding quantile of observation \code{y}
  together with the median of the predictive distribution.
}
\keyword{regression}
