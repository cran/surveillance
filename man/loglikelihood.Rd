\name{loglikelihood}
\alias{loglikelihood}
\encoding{latin1}
\title{Calculation of the loglikelihood needed in algo.hhh}
\description{
        Calculates the loglikelihood according to the model specified in \code{designRes}.
}
\usage{
loglikelihood(theta, designRes)
}
\arguments{
        \item{theta}{vector of parameters 
              \deqn{\theta = (\alpha_1,\ldots,\alpha_m,
                              \bold{\lambda}, 
                              \bold{\phi}, \bold{\beta}, 
										\bold{\gamma}_1, \ldots,  \bold{\gamma}_m, 
										\bold{\psi}),}
				  where \eqn{\bold{\lambda}=(\lambda_1,\ldots,\lambda_m)}, 
				  \eqn{\bold{\phi}=(\phi_1,\ldots,\phi_m)}, 
				  \eqn{\bold{\beta}=(\beta_1,\ldots,\beta_m)}, 
				  \eqn{\bold{\gamma_1}=(\gamma_{11},\ldots,\gamma_{1,2S_1})}{\bold{\gamma_1}=(\gamma_11,\ldots,\gamma_(1,2S_1))}, 
				  \eqn{\bold{\gamma_m}=(\gamma_{m1},\ldots,\gamma_{m,2S_m})}{\bold{\gamma_m}=(\gamma_m1,\ldots,\gamma_(m,2S_m))}, 
				  \eqn{\bold{\psi}=(\psi_1,\ldots,\psi_m)}.
				  
              If the model specifies less parameters, those components are omitted.} 
        \item{designRes}{Result of a call to \code{make.design} }
}
\value{
    Returns the loglikelihood
}
\seealso{\code{\link{meanResponse}}}
\author{M. Paul, L. Held}
\keyword{models}


% meanResponse
