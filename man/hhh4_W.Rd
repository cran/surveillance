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
  \eqn{w_{ji} = o_{ji}^{-d}}{w_ji = o_ji^-d}
  (if \eqn{o_{ji} > 0}{o_ji > 0}, otherwise \eqn{w_{ji} = 0}{w_ji = 0}),
  where \eqn{o_{ji}}{o_ji} (\eqn{=o_{ij}}{=o_ij}) is the adjacency order
  between regions \eqn{i} and \eqn{j},
  and the decay parameter \eqn{d} is to be estimated.
  In the nonparametric formulation, unconstrained log-weights will be
  estimated for each of the adjacency orders \code{2:maxlag} (the
  first-order weight is fixed to 1 for identifiability).
  Both weight functions can be modified to include a 0-distance weight,
  which enables \code{hhh4} models without a separate autoregressive component.
}

\usage{
W_powerlaw(maxlag, normalize = TRUE, log = FALSE,
           initial = if (log) 0 else 1, from0 = FALSE)

W_np(maxlag, truncate = TRUE, normalize = TRUE,
     initial = log(zetaweights(2:(maxlag+from0))),
     from0 = FALSE, to0 = truncate)
}

\arguments{
  \item{maxlag}{a single integer specifying a limiting order of
    adjacency. If spatial dependence is not to be truncated at some
    high order, \code{maxlag} should be set to the maximum adjacency
    order in the network of regions. The smallest possible value for
    \code{maxlag} is 2 if \code{from0=FALSE} and 1 otherwise.}
  \item{truncate,to0}{\code{W_np} represents order-specific log-weights up to
    order \code{maxlag}. Higher orders are by default (\code{truncate=TRUE})
    assumed to have zero weight (similar to \code{W_powerlaw}).
    Alternatively, \code{truncate=FALSE} requests that the weight at
    order \code{maxlag} should be carried forward to higher orders.
    \code{truncate} has previously been called \code{to0} (deprecated).}
  \item{normalize}{logical indicating if the weights should be normalized
    such that the rows of the weight matrix sum to 1 (default).
    Note that normalization does not work with islands, i.e., regions
    without neighbours.}
  \item{log}{logical indicating if the decay parameter \eqn{d} should be
    estimated on the log-scale to ensure positivity.}
  \item{initial}{initial value of the parameter vector.}
  \item{from0}{logical indicating if these parametric weights should
    include the 0-distance (autoregressive) case. In the default setting
    (\code{from0 = FALSE}), adjacency order 0 has zero weight, which is
    suitable for \code{hhh4} models with a separate autoregressive
    component. With \code{from0 = TRUE} (Meyer and Held, 2017), the
    power law is based on \eqn{(o_{ji} + 1)}{(o_ji + 1)}, and
    nonparametric weights are estimated for adjacency orders
    \code{1:maxlag}, respectively, where the 0-distance weight is
    \eqn{w_{jj} = 1}{w_jj = 1} (without normalization). Note that
    the corresponding \code{hhh4} model should then exclude a separate
    autoregressive component (\code{control$ar$f = ~ -1}).}
}

\value{
  a list which can be passed as a specification of parametric
  neighbourhood weights in the \code{control$ne$weights} argument of
  \code{\link{hhh4}}.
}

\details{
  \code{hhh4} will take adjacency orders from the \code{neighbourhood}
  slot of the \code{"sts"} object, so these must be prepared before
  fitting a model with parametric neighbourhood weights. The function
  \code{\link{nbOrder}} can be used to derive adjacency orders from a
  binary adjacency matrix.
}

\references{
  Meyer, S. and Held, L. (2014):
  Power-law models for infectious disease spread.
  \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.
  \doi{10.1214/14-AOAS743}

  Meyer, S. and Held, L. (2017):
  Incorporating social contact data in spatio-temporal models for
  infectious disease spread.
  \emph{Biostatistics}, \bold{18} (2), 338-351.
  \doi{10.1093/biostatistics/kxw051}
}

\author{
  Sebastian Meyer
}

\seealso{  
  \code{\link{nbOrder}} to determine adjacency orders from a binary
  adjacency matrix.
  
  \code{\link{getNEweights}} and \code{\link{coefW}} to extract the
  estimated neighbourhood weight matrix and coefficients from an
  \code{hhh4} model.
}

\examples{
data("measlesWeserEms")

## data contains adjaceny orders as required for parametric weights
plot(measlesWeserEms, type = observed ~ unit, labels = TRUE)
neighbourhood(measlesWeserEms)[1:6,1:6]
max(neighbourhood(measlesWeserEms))  # max order is 5

## fit a power-law decay of spatial interaction
## in a hhh4 model with seasonality and random intercepts in the endemic part
measlesModel <- list(
    ar = list(f = ~ 1),
    ne = list(f = ~ 1, weights = W_powerlaw(maxlag=5)),
    end = list(f = addSeason2formula(~-1 + ri(), S=1, period=52)),
    family = "NegBin1")

## fit the model
set.seed(1)  # random intercepts are initialized randomly
measlesFit <- hhh4(measlesWeserEms, measlesModel)
summary(measlesFit)  # "neweights.d" is the decay parameter d
coefW(measlesFit)

## plot the spatio-temporal weights o_ji^-d / sum_k o_jk^-d
## as a function of adjacency order
plot(measlesFit, type = "neweights", xlab = "adjacency order")
## normalization => same distance does not necessarily mean same weight.
## to extract the whole weight matrix W: getNEweights(measlesFit)

## visualize contributions of the three model components
## to the overall number of infections (aggregated over all districts)
plot(measlesFit, total = TRUE)
## little contribution from neighbouring districts

## simpler model with autoregressive effects captured by the ne component
measlesModel2 <- list(
    ne = list(f = ~ 1, weights = W_powerlaw(maxlag=5, from0=TRUE)),
    end = list(f = addSeason2formula(~-1 + ri(), S=1, period=52)),
    family = "NegBin1")
measlesFit2 <- hhh4(measlesWeserEms, measlesModel2)
## omitting the separate AR component simplifies model extensions/selection
## and interpretation of covariate effects (only two predictors left)

plot(measlesFit2, type = "neweights", exclude = NULL, xlab = "adjacency order")
## strong decay, again mostly within-district transmission
## (one could also try a purely autoregressive model)
plot(measlesFit2, total = TRUE,
     legend.args = list(legend = c("epidemic", "endemic")))
## almost the same RMSE as with separate AR and NE effects
sqrt(mean(residuals(measlesFit, "response")^2))
sqrt(mean(residuals(measlesFit2, "response")^2))
}

\keyword{spatial}
\keyword{models}
\keyword{utilities}
