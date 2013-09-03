\name{powerlaw}
\alias{powerlaw}
\title{
  Power-Law Neighbourhood Weights for \code{hhh4} Models
}

\description{
  Set up a power-law weight structure for the neighbourhood
  component of \code{\link{hhh4}} models, i.e., without normalization
  we have \eqn{w_{ji} = o_{ji}^{-d}}{w_ji = o_ji^-d}, where
  \eqn{o_{ji}}{o_ji} is the order of neighbourhood between regions
  \eqn{i} and \eqn{j}, see \code{\link{nbOrder}}.
}

\usage{
powerlaw(maxlag, normalize = TRUE, log = FALSE, initial = if (log) 0 else 1)
}

\arguments{
  \item{maxlag}{single integer specifying up to which order of
    neighbourhood the weights should be positive; they will be set to 0
    for higher orders. If spatial dependence is not to be truncated at
    some high order, \code{maxlag} should be set to the maximum
    neighbourhood order in the network of regions.}
  \item{normalize}{logical indicating if the weights should be normalized
    such that the rows of the weight matrix sum to 1.}
  \item{log}{logical indicating if the decay parameter should be
    estimated on the log-scale to ensure positivity.}
  \item{initial}{scalar initial value for the (log-)decay parameter.}
}

\value{
  a list which can be passed as a specification of parametric
  neighbourhood weights in \code{\link{hhh4}}.
}

\references{
  Meyer, S. and Held, L. (2013):
  Modelling power-law spread of infectious diseases.
  \emph{arXiv:1308.5115}\cr
  arXiv-Link: \url{http://arxiv.org/abs/1308.5115}
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{nbOrder}} to determine the matrix of neighbourhood orders
  from a binary adjacency matrix.
  
  \code{\link{siaf.powerlaw}} for modelling power-law distance decay in
  \code{\link{twinstim}} space-time point process models.
}

\examples{
if (requireNamespace("spdep")) {
  data("measles.weser")
  measles <- disProg2sts(measles.weser)

  ## parametric powerlaw weights require neighbourhood orders in the data
  neighbourhood(measles) <- nbOrder(neighbourhood(measles), maxlag=Inf)

  ## a simple hhh4 model with power-law decay of spatial interaction
  m <- list(ar = list(f = ~ 1),
            ne = list(f = ~ 1, weights = powerlaw(maxlag=4, log=FALSE)),
            end = list(f = addSeason2formula(~-1 + ri(), S=1, period=52),
                       offset = population(measles)),
            family = "NegBin1", verbose=TRUE)

  ## fit the model
  set.seed(1)  # random intercepts are initialized randomly
  fit <- hhh4(measles, m)
  summary(fit)  # "neweights.d" is the decay parameter
}
}

\keyword{spatial}
\keyword{models}
\keyword{utilities}
