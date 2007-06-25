\name{algo.hhh}
\alias{algo.hhh}
\alias{print.ah}
\title{Model fit based on the Held, Höhle, Hofmann paper}
\description{
        Fits a Poisson/Negative Binomial model with mean \eqn{\mu_{it}} 
        (as described in Held/Höhle/Hofmann, 2005) to a multivariate 
        time series of counts.  
}
\usage{
  algo.hhh(disProgObj, control=list(linear=FALSE, nseason=0, period=52,   
      neighbours=FALSE, negbin=FALSE, lambda=TRUE), thetastart=NULL, 
      verbose=TRUE)
}
\arguments{
\item{disProgObj}{object of class \code{disProg}} 
\item{control}{Control object:
    \itemize{
    \item{linear}{if \code{true} a linear trend is included.}
    \item{nseason}{number of Fourier frequencies to be used}
    \item{period}{period length (defaults to \code{disProgObj$freq}, i.e. 52 for weekly data)}
    \item{neighbours}{if \code{true} a autoregressive parameter for 
        adjacent units is included.}
    \item{negbin}{if \code{true} negative binomial rather than poisson is used.}
    \item{lambda}{if \code{true} a autoregressive parameter is included.}
    }}
\item{thetastart}{vector with starting values for all parameters specified
    in the control object (for \code{optim}).  } 
\item{verbose}{if \code{true} information about convergence is printed}
}
\value{\item{list}{
    Returns a list with estimated parameters, standard errors, 
    covariance matrix and the value of the loglikelihood as elements. 
  }
}
\details{Note that for the time being this function is not a surveillance algorithm, but only a modelling approach as described in the Held et. al (2005) paper.}
\seealso{\code{\link{meanResponse}}}
\author{M. Paul, L. Held, M. Höhle}
\examples{

    # univariate time series: salmonella agona cases
    data(salmonella.agona)
    salmonella <- create.disProg(week=1:length(salmonella.agona),
                 observed=salmonella.agona$observed, state=salmonella.agona$state)
    
    model1 <- list(lambda=TRUE, linear=TRUE, nseason=1, negbin=TRUE)
    algo.hhh(salmonella, control=model1)
    
    # multivariate time series: measels cases
    data(measels.weser)
    
    # same model as above
    algo.hhh(measels.weser, control=model1)
    
    # different starting values for (lambda, beta, gamma, delta, psi)
    startValues <- c(0.1, rep(0, 3), 1)
    algo.hhh(measels.weser, control=model1, thetastart=startValues)
    
    # include autoregressive parameter phi for adjacent "Kreise"
    model2 <- list(lambda=TRUE, neighbours=TRUE, linear=TRUE, nseason=1,
                  negbin=TRUE)
    algo.hhh(measels.weser, control=model2)
                  
}
\keyword{}
\source{Held, L., Höhle, M., Hofmann, M. (2005). A statistical framework for the 
    analysis of multivariate infectious disease surveillance counts. Statistical 
    Modelling, 5, p. 187-199.}


% predict.ah
