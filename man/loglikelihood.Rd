\name{loglikelihood}
\alias{loglikelihood}
\title{Calculation of the loglikelihood needed in algo.hhh}
\description{
        Calculates the loglikelihood according to the model specified in \code{designRes}.
}
\usage{
  loglikelihood(theta, designRes)
}
\arguments{
        \item{theta}{vector of parameters 
              \eqn{\theta = (\lambda, \phi, \beta, \gamma_1, \delta_1, 
              \gamma_2, \delta_2, \ldots, \psi, \alpha_1, \alpha_2, \ldots)}.  
              If the model specifies less parameters, those components are omitted.}
        \item{designRes}{Result of a call to \code{make.design}}
}
\value{
    Returns the loglikelihood
}
\details{}
\seealso{}
\author{M. Paul, L. Held}
\examples{}
\keyword{}
\source{}



% meanResponse
