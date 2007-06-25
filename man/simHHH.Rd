\name{simHHH}
\alias{simHHH}
\alias{simHHH.default}
\alias{simHHH.ah}
\title{Simulates data based on the model proposed by Held et. al (2005)}
\description{
   Simulates a multivariate time series of counts based on the
   Poisson/Negative Binomial model as described in Held et al. (2005).
}
\usage{
\method{simHHH}{default}(model=NULL, control = list(coefs = list(alpha, gamma = 0, delta = 0,
       lambda = 0, phi = NULL, psi = NULL, period = 52),
       neighbourhood = NULL, population = NULL, start = NULL),
%       state = NULL, K = 0),
       length)

\method{simHHH}{ah}(model, control = model$control, length)
}
\arguments{
        \item{control}{list with
             \describe{
             \item{coefs}{list with the following parameters of the model -
                  if not specified, those parameters are omitted
                   \describe{
                      \item{alpha}{vector of length \code{m} with intercepts
                                   for \code{m} units or geographic areas respectively}
                      \item{gamma}{vector with parameters for the "sine" part of \eqn{\nu_{i,t}} }
                      \item{delta}{vector with parameters for the "cosine" part of \eqn{\nu_{i,t}} }
                      \item{lambda}{autoregressive parameter}
                      \item{phi}{autoregressive parameter for adjacent units}
                      \item{psi}{overdispersion parameter of the negative binomial
                            model; \code{NULL} corresponds to a Poisson model}
                      \item{period}{period of the seasonal component, defaults
                            to 52 for weekly data}
                      }  }
             \item{neighbourhood}{neighbourhood matrix of size
                  \eqn{m \times m} with element 1 if two units are adjacent;
                  the default \code{NULL} assumes that there are no neighbours}
             \item{population}{matrix with population proportions;
                     the default \code{NULL} sets \eqn{n_{i,t}=1} }
             \item{start}{if \code{NULL}, the means of the endemic part in the
                   \code{m} units is used  as initial values \eqn{y_{i,0}} }
             }}
        \item{model}{Result of a model fit with
                     \code{\link{algo.hhh}}, the estimated parameters
                      are used to simulate data}
        \item{length}{number of time points to simulate  }
}
\value{Returns a list with elements
       \item{data}{\code{disProgObj} of simulated data }
       \item{mean}{matrix with mean \eqn{\mu_{i,t}} that was used to
                   simulate the data}
       \item{endemic}{matrix with only the endemic part \eqn{\nu_{i,t}}
        }
       \item{coefs}{list with parameters of the model}
}
\details{
        Simulates data from a Poisson or a Negative Binomial model
        with mean
        \deqn{\mu_{it} = \lambda y_{i,t-1} + \phi \sum_{j \sim i} y_{j,t-1} + n_{it} \nu_{it}}
        where
        \deqn{\log \nu_{it} = \alpha_i
        + \sum_{s=1}^{S}(\gamma_s sin(\omega_s t) + \delta_s cos(\omega_s t))}
        \eqn{\omega_s = 2s\pi/}\code{period} are Fourier frequencies
        and \eqn{n_{it}} are possibly standardized population sizes.
}
\note{The model does not contain a linear trend.}
%\seealso{\code{\link{algo.hhh}}, \code{\link{algo.hhh.grid}}}
\source{Held, L., Höhle, M., Hofmann, M. (2005). A statistical framework for the
    analysis of multivariate infectious disease surveillance counts. Statistical
    Modelling, 5, p. 187-199.
    }
\keyword{}

