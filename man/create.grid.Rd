\name{create.grid}
\alias{create.grid}
\title{Computes a matrix of initial values}
\description{
       For a given list of parameters specified as \code{param = c(lower,upper,length)}
       \code{create.grid} creates a grid of initial values for \code{algo.hhh.grid}.
       The resulting matrix contains all combinations of the supplied parameters 
       which each are a sequence of length \code{length} from \code{lower} to 
       \code{upper}.  
       Note that \eqn{\lambda, \phi, \psi} must be positive.  
       For \eqn{\gamma, \delta} only one sequence of starting
       values is considered to create the grid, i.e. intital values are the same
       for each one of the seasonal parameters. \code{nseason} is the number of 
       Fourier frequencies to be included.       
}
\usage{
  create.grid(params = list(lambda=c(0.1,0.9,5), phi=c(0.1,0.9,5), 
      psi=c(0.3,12,10), beta=c(-0.5,0.5,3), gammaDelta=c(-0.5,0.5,3), 
      nseason=1))
}
\arguments{
        \item{params}{list of parameters: \code{param=c(lower,upper,length)}
            \itemize{
            \item{lambda}{autoregressive parameter \eqn{\lambda}{lambda}. }
            \item{phi}{autoregressive parameter for adjacent units \eqn{\phi}{phi}. }
            \item{psi}{overdispersion parameter for Negative Binomial model \eqn{\psi}{psi}. }
            \item{beta}{linear trend \eqn{\beta}{beta}.}
            \item{gammaDelta}{seasonal parameters \eqn{\gamma, \delta}. }
            \item{nseason}{number of Fourier frequencies \eqn{S} to be used  
                (if not specified, \code{nseason} defaults to 1)}
            }
          }
}
\value{\item{matrix}{matrix with \code{gridSize} starting values as rows}
       
}
\details{}
\seealso{\code{\link{algo.hhh.grid}}}
\author{M. Paul}
\examples{
    
    # consider the model specified in a control object for algo.hhh.grid
    cntrl1 <- list(lambda=TRUE, neighbours=TRUE, linear=TRUE, nseason=1)
    cntrl2 <- list(lambda=TRUE, negbin=TRUE, linear=TRUE)
    
    # create a grid of initial values for respective parameters
    grid1 <- create.grid(params = list(lambda=c(0.1,0.9,5), phi=c(0.1,0.9,5), 
                         beta=c(-1,1,3), gammaDelta=c(-1,1,3), nseason=1))
    grid2 <- create.grid(params = list(lambda=c(0.1,0.9,5), psi=c(0.3,12,10),
                         beta=c(0,0,1)))

}
\keyword{}
\source{}


% loglikelihood
