\name{meanResponse}
\alias{meanResponse}
\title{Calculate mean response needed in algo.hhh}
\description{
        Calculates the mean response for the model specified in designRes 
        according to equations (1.2) and (1.1) in Held et al., 2005
        for univariate time series and equations (3.3) and (3.2) for multivariate
        time series. See details.
}
\usage{
  meanResponse(theta, designRes)
}
\arguments{
        \item{theta}{vector of parameters 
              \deqn{\theta = (\lambda, \phi, \beta, \gamma_1, \delta_1, 
              \gamma_2, \delta_2, \ldots, \psi, \alpha_1, \alpha_2, \ldots).}
              If the model specifies less parameters, those components are omitted.} 
        \item{designRes}{Result of a call to \code{make.design}}
}
\value{
    Returns a matrix of dimension \eqn{n \times m} with the calculated mean 
    response for each time point and unit, where \eqn{n} is the number of time 
    points and \eqn{m} is the number of units.
}
\details{        
        Calculates the mean response for a Poisson or a Negative Binomial model
        with mean
        \deqn{\mu_t = \nu_t + \lambda y_{t-1}}
        where 
        \deqn{log \nu_t = \alpha + \beta t 
        + \sum_{s=1}^{S}(\gamma_s sin(\omega_s t) + \delta_s cos(\omega_s t))}
        and Fourier frequencies \eqn{\omega_s = 2s\pi/period} for a univariate time series.  
        For multivariate time series the mean structure is
        \deqn{\mu_{it} = \lambda y_{i,t-1} + \phi \sum_{j \sim i} y_{j,t-1} + n_{it} \nu_{it}}
        where 
        \deqn{log \nu_{it} = \alpha_i + \beta t 
        + \sum_{s=1}^{S}(\gamma_s sin(\omega_s t) + \delta_s cos(\omega_s t))}
        and \eqn{n_{it}} are standardized population counts.
        }
\seealso{}
\author{M. Paul, L. Held}
\examples{}
\keyword{}
\source{Held, L., Höhle, M., Hofmann, M. (2005). A statistical framework for the 
    analysis of multivariate infectious disease surveillance counts. Statistical 
    Modelling, 5, p. 187-199.}


% make.design
