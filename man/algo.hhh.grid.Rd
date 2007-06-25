\name{algo.hhh.grid}
\alias{algo.hhh.grid}
\alias{print.ahg}
\title{Function to try multiple starting values}}
\description{
        Tries multiple starting values in \code{algo.hhh}. 
        Starting values are provided in a matrix with \code{gridSize} rows, the
        grid search is conducted until either all starting values are used or
        a time limit \code{maxTime} is exceeded.
        The result with the highest likelihood is returned. 
}
\usage{
  algo.hhh.grid(disProgObj, control=list(linear=FALSE, nseason=0, period=52,   
      neighbours=FALSE, negbin=FALSE, lambda=TRUE), thetastartMatrix,
      maxTime=1800, verbose=FALSE)
}
\arguments{
        \item{disProgObj}{object of class \code{disProg}} 
        \item{control}{Control object:
             \itemize{
             \item{linear}{if \code{true} a linear trend is included.}
             \item{nseason}{number of Fourier frequencies to be used}
             \item{period}{period length (defaults to \code{disProgObj$freq}, i.e. 52 for weekly data)}
             \item{neighbours}{if \code{true} a autoregressive parameter for 
                 adjacent areas is included.}
             \item{negbin}{if \code{true} negative binomial rather than poisson is used.}
             \item{lambda}{if \code{true} a autoregressive parameter is included.}
             }}
        \item{thetastartMatrix}{matrix with initial values for all parameters specified
             in the control object as rows.}
        \item{verbose}{if \code{true} progress information is printed} 
        \item{maxTime}{maximum of time (in seconds) to elapse until algorithm stopps.}
}
\value{\item{list}{
        \itemize{
             \item{best}{result of a call to \code{algo.hhh} with highest likelihood } 
             \item{allLoglik}{values of loglikelihood for all starting values used}
             \item{gridSize}{number of different starting values in thetastartMatrix}
             \item{gridUsed}{number of used starting values}
             \item{convergence}{if \code{false} \code{algo.hhh} did not converge
                 for all (used) starting values}
        }} 
}
\details{}
\seealso{\code{\link{meanResponse}},\code{\link{create.grid}},\code{\link{algo.hhh}}}
\author{M. Paul, L. Held}
\examples{
    
    # monthly counts of menigococcal infections in France
    data(meningo.age)
                  
    # specify model for algo.hhh.grid
    model1 <- list(lambda=TRUE, linear=TRUE)
    # create grid of inital values
    grid1 <- create.grid(params = list(lambda=c(0.1,0.9,5), beta=c(-1,1,3)))
    # try multiple starting values, print progress information
    algo.hhh.grid(meningo.age, control=model1, thetastartMatrix=grid1,
                  verbose=TRUE)
    
    # specify model
    model2 <- list(lambda=TRUE, neighbours=TRUE, linear=TRUE, negbin=TRUE,
                  nseason=1, period=12)
    grid2 <- create.grid(params = list(lambda=c(0.1,0.9,5), phi=c(0.1,0.9,5), 
                  beta=c(0,0.5,2), psi=c(0.3,12,10), gammaDelta=c(-0.5,0.5,3),
                  nseason=1))
    # run algo.hhh.grid, search time is limited to 60 sec
    algo.hhh.grid(meningo.age, control=model2, thetastartMatrix=grid2,
                  maxTime=60)

}
\keyword{}
\source{Held, L., Höhle, M., Hofmann, M. (2005). A statistical framework for the 
    analysis of multivariate infectious disease surveillance counts. Statistical 
    Modelling, 5, p. 187-199.}


% create.grid
