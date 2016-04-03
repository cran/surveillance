\name{algo.hhh.grid}
\alias{algo.hhh.grid}
\alias{print.ahg}
\alias{coef.ahg}
\encoding{latin1}
\title{Fit a Classical HHH Model (DEPRECATED) with Varying Start Values}
\description{
        \code{algo.hhh.grid} tries multiple starting values in \code{algo.hhh}. 
        Starting values are provided in a matrix with \code{gridSize} rows
        (usually created by \code{\link{create.grid}}).
        The grid search is conducted until either all starting values are used or
        a time limit \code{maxTime} is exceeded.
        The result with the highest likelihood is returned.
        
        Note that the \code{algo.hhh} implementation of HHH models is
        \strong{deprecated} and superseded by the function \code{\link{hhh4}}.
}
\usage{
algo.hhh.grid(disProgObj, control=list(lambda=TRUE, neighbours=FALSE, 
   linear=FALSE, nseason=0,
   negbin=c("none", "single", "multiple"), 
   proportion=c("none", "single", "multiple"),lag.range=NULL), 
   thetastartMatrix, maxTime=1800, verbose=FALSE)
}
\arguments{
\item{disProgObj}{object of class \code{disProg}}
\item{control}{control object:
      \describe{
      \item{\code{lambda}}{If \code{TRUE} an autoregressive 
            parameter \eqn{\lambda} is included, if \code{lambda} is a vector of logicals,
            unit-specific parameters \eqn{\lambda_i} are included. By default, observations \eqn{y_{t-lag}}{y_t-lag}
            at the previous time points, i.e. \eqn{lag=1}, are used for the autoregression. Other lags can be
            used by specifying \code{lambda} as a vector of integers, see Examples and \code{\link{algo.hhh}} for details.}
      \item{\code{neighbours}}{If \code{TRUE} an autoregressive parameter for 
            adjacent units \eqn{\phi} is included, if \code{neighbours} is a vector of logicals,
            unit-specific parameters \eqn{\phi_i} are included. By default, observations \eqn{y_{t-lag}}{y_t-lag}
            at the previous time points, i.e. \eqn{lag=1}, are used for the autoregression. Other lags can be
            used by specifying \code{neighbours} as a vector of integers.}
      \item{\code{linear}}{a \code{logical} (or a vector of logicals) indicating wether a linear 
            trend \eqn{\beta} (or a linear trend \eqn{\beta_i} for each unit)
            is included}
      \item{\code{nseason}}{integer number of Fourier frequencies; if \code{nseason} is a vector
            of integers, each unit \eqn{i} gets its own seasonal parameters }
      \item{\code{negbin}}{if \code{"single"} negative binomial rather than poisson is used,
            if \code{"multiple"} unit-specific overdispersion parameters are used.}
      \item{\code{proportion}}{see details in \code{\link{algo.hhh}} }
      \item{\code{lag.range}}{determines which observations are used to fit the model}
      }}        
      \item{thetastartMatrix}{matrix with initial values for all parameters specified
      in the control object as rows.}
\item{verbose}{if \code{true} progress information is printed} 
\item{maxTime}{maximum of time (in seconds) to elapse until algorithm stopps.}
}

\value{Returns an object of class \code{ahg} with elements 
    \item{best}{result of a call to \code{algo.hhh} with highest likelihood } 
    \item{allLoglik}{values of loglikelihood for all starting values used}
    \item{gridSize}{number of different starting values in thetastartMatrix}
    \item{gridUsed}{number of used starting values}
    \item{time}{elapsed time}
    \item{convergence}{if \code{false} \code{algo.hhh} did not converge
        for all (used) starting values}
}

\seealso{\code{\link{algo.hhh}}, \code{\link{create.grid}}
}
\author{M. Paul, L. Held}
\examples{
## monthly counts of menigococcal infections in France
data(meningo.age)
            
# specify model for algo.hhh.grid
model1 <- list(lambda=TRUE)

# create grid of inital values
grid1 <- create.grid(meningo.age, model1, 
                     params = list(epidemic=c(0.1,0.9,5)))
                     
# try multiple starting values, print progress information
algo.hhh.grid(meningo.age, control=model1, thetastartMatrix=grid1,
               verbose=TRUE)


## more sophisticated models with a much longer runtime follow
\dontrun{
# specify model
model2 <- list(lambda=TRUE, neighbours=TRUE, negbin="single",
               nseason=1)
grid2 <- create.grid(meningo.age, model2, 
                     params = list(epidemic=c(0.1,0.9,3), 
                                   endemic=c(-0.5,0.5,3), 
                                   negbin = c(0.3, 12, 10)))
                                 
# run algo.hhh.grid, search time is limited to 30 sec
algo.hhh.grid(meningo.age, control=model2, thetastartMatrix=grid2,
              maxTime=30)
              
## weekly counts of influenza and meningococcal infections in Germany, 2001-2006
data(influMen)

# specify model with two autoregressive parameters lambda_i, overdispersion
# parameters psi_i, an autoregressive parameter phi for meningococcal infections
# (i.e. nu_flu,t = lambda_flu * y_flu,t-1  
#  and  nu_men,t = lambda_men * y_men,t-1 + phi_men*y_flu,t-1 )
# and S=(3,1) Fourier frequencies
model <- list(lambda=c(TRUE,TRUE), neighbours=c(FALSE,TRUE),
              linear=FALSE, nseason=c(3,1),negbin="multiple")
              
# create grid of initial values
grid <- create.grid(influMen,model, list(epidemic=c(.1,.9,3),
              endemic=c(-.5,.5,3), negbin=c(.3,15,10)))
# run algo.hhh.grid, search time is limited to 30 sec
algo.hhh.grid(influMen, control=model, thetastartMatrix=grid, maxTime=30)

# now meningococcal infections in the same week should enter as covariates
# (i.e. nu_flu,t = lambda_flu * y_flu,t-1  
#  and  nu_men,t = lambda_men * y_men,t-1 + phi_men*y_flu,t )
model2 <- list(lambda=c(1,1), neighbours=c(NA,0),
              linear=FALSE,nseason=c(3,1),negbin="multiple")
              
algo.hhh.grid(influMen, control=model2, thetastartMatrix=grid, maxTime=30)
}
}
\keyword{ts}
\keyword{regression}
\references{
  Held, L., \enc{Höhle}{Hoehle}, M., Hofmann, M. (2005) A statistical framework 
    for the analysis of multivariate infectious disease surveillance counts, 
    Statistical Modelling, \bold{5}, 187--199.
    
  Paul, M., Held, L. and Toschke, A. M. (2008) Multivariate modelling of 
    infectious disease surveillance data, Statistics in Medicine, \bold{27}, 
    6250--6267.    
}
