\name{algo.farrington.threshold}
\alias{algo.farrington.threshold}
\encoding{latin1}

\title{Threshold computations using a two sided confidence interval}
\description{
Depending on the current transformation \eqn{h(y)= \{y, \sqrt{y}, y^{2/3}\}},

  \deqn{V(h(y_0)-h(\mu_0))=V(h(y_0))+V(h(\mu_0))}

  is used to compute a prediction interval. The prediction variance
  consists of a component due to the variance of having a single
  observation and a prediction variance.  }

\usage{
algo.farrington.threshold(pred,phi,alpha=0.01,skewness.transform="none")
}
\arguments{
\item{pred}{A GLM prediction object}
\item{phi}{Current overdispersion (superflous?)}
\item{alpha}{Quantile level in Gaussian based CI, i.e. an \eqn{(1-\alpha)\%}
    confidence interval is computed. }
\item{skewness.transform}{Skewness correction, i.e. one of
    \code{"none"}, \code{"sqrt"}, or \code{"2/3"}.}
}
\value{\item{vector}{
        Vector of length 2 with lower and upper bounds of an \eqn{(1-\alpha)\%} confidence interval.
        }
}
%\keyword{}


