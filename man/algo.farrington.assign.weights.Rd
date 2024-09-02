\name{algo.farrington.assign.weights}
\alias{algo.farrington.assign.weights}
\title{Assign weights to base counts}
\description{
  Weights are assigned according to the Anscombe residuals
}
\usage{
algo.farrington.assign.weights(s, weightsThreshold=1)
}
\arguments{
  \item{s}{Vector of standardized Anscombe residuals}
  \item{weightsThreshold}{A scalar indicating when observations are seen
    as outlier. In the original Farrington proposal the value was 1
    (default value), in the improved version this value is suggested to be 2.58.}
}
\value{Weights according to the residuals}
\seealso{\code{\link{anscombe.residuals}}}
\keyword{regression}
