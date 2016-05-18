\name{hhh4_W}
\alias{W_powerlaw}
\alias{W_np}
\title{
  Power-Law and Nonparametric Neighbourhood Weights for \code{hhh4}-Models
}

\description{
  Set up power-law or nonparametric weights for the neighbourhood
  component of \code{\link{hhh4}}-models as proposed by Meyer and Held (2014).
  Without normalization, power-law weights are
  \eqn{w_{ji} = o_{ji}^{-d}}{w_ji = o_ji^-d}, where \eqn{o_{ji}}{o_ji}
  is the order of neighbourhood between regions \eqn{i} and \eqn{j},
  see \code{\link{nbOrder}}, and \eqn{d} is to be estimated.
  In the nonparametric formulation, \code{maxlag-1} order-specific
  log-weights are to be estimated (the first-order weight is always
  fixed to 1 for identifiability).

}

\usage{
W_powerlaw(maxlag, normalize = TRUE,
           log = FALSE, initial = if (log) 0 else 1)

W_np(maxlag, to0 = TRUE, normalize = TRUE,
     initial = log(zetaweights(2:maxlag)))
}

\arguments{
  \item{maxlag}{a single integer specifying a limiting order of
    neighbourhood. If spatial dependence is not to be truncated at some
    high order, \code{maxlag} should be set to the maximum neighbourhood
    order in the network of regions.}
  \item{to0}{\code{W_np} represents order-specific log-weights up to
    order \code{maxlag}. Higher orders are by default (\code{to0=TRUE})
    assumed to have 0 weight as for \code{W_powerlaw}. Alternatively,
    \code{to0=FALSE} requests that the weight at order \code{maxlag}
    should be carried forward to higher orders.}
  \item{normalize}{logical indicating if the weights should be normalized
    such that the rows of the weight matrix sum to 1 (default).
    Note that normalization does not work with islands, i.e., regions
    without neighbours.}
  \item{log}{logical indicating if the decay parameter \eqn{d} should be
    estimated on the log-scale to ensure positivity.}
  \item{initial}{initial value of the parameter vector.}
}

\value{
  a list which can be passed as a specification of parametric
  neighbourhood weights in the \code{control$ne$weights} argument of
  \code{\link{hhh4}}.
}

\references{
  Meyer, S. and Held, L. (2014):
  Power-law models for infectious disease spread.
  \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.
  \doi{10.1214/14-AOAS743}
}

\author{
  Sebastian Meyer
}

\seealso{  
  \code{\link{nbOrder}} to determine the matrix of neighbourhood orders
  from a binary adjacency matrix.
  
  \code{\link{siaf.powerlaw}}, and \code{\link{siaf.step}} for modelling
  distance decay as power law or step function in 
  \code{\link{twinstim}} space-time point process models.
}

\examples{
data("measlesWeserEms")

## data contains neighbourhood orders as required for parametric weights
neighbourhood(measlesWeserEms)[1:6,1:6]
max(neighbourhood(measlesWeserEms))  # max order is 5

## fit a power-law decay of spatial interaction
## in a hhh4 model with seasonality and random intercepts in the endemic part
measlesModel <- list(
    ar = list(f = ~ 1),
    ne = list(f = ~ 1, weights = W_powerlaw(maxlag=5, normalize=TRUE, log=FALSE)),
    end = list(f = addSeason2formula(~-1 + ri(), S=1, period=52),
               offset = population(measlesWeserEms)),
    family = "NegBin1")

## fit the model
set.seed(1)  # random intercepts are initialized randomly
measlesFit <- hhh4(measlesWeserEms, measlesModel)
summary(measlesFit)  # "neweights.d" is the decay parameter d

## plot the spatio-temporal weights o_ji^-d / sum_k o_jk^-d
## as a function of neighbourhood order
plot(measlesFit, type="neweights")
## Due to normalization, same distance does not necessarily mean same weight.
## There is no evidence for a power law of spatial interaction in this
## small observation region with only 17 districts.
## A possible simpler model is first-order dependence, i.e., using
## 'weights = neighbourhood(measlesWeserEms) == 1' in the 'ne' component.
}

\keyword{spatial}
\keyword{models}
\keyword{utilities}
