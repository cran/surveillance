\name{algo.hhh}
\alias{algo.hhh}
\alias{print.ah}
\alias{coef.ah}
\encoding{latin1}
\title{Model fit based on the Held, Höhle, Hofman paper} 
\description{
        Fits a Poisson/negative binomial model with mean \eqn{\mu_{it}}{\mu_it} 
        (as described in Held/Höhle/Hofmann, 2005) to a multivariate 
        time series of counts.  
}
\usage{
algo.hhh(disProgObj, control=list(lambda=TRUE, neighbours=FALSE, 
   linear=FALSE, nseason = 0,
   negbin=c("none", "single", "multiple"), 
   proportion=c("none", "single", "multiple")), 
   thetastart=NULL, verbose=TRUE)
}
\arguments{
        \item{disProgObj}{Object of class \code{disProg}}
        \item{control}{Control object:
              \describe{
              \item{\code{lambda}}{if \code{TRUE} an autoregressive 
                   parameter \eqn{\lambda} is included, if \code{lambda} is a vector of logicals,
                   unit-specific parameters \eqn{\lambda_i} are included}
              \item{\code{neighbours}}{if \code{TRUE} an autoregressive parameter for 
                    adjacent units \eqn{\phi} is included, if \code{neighbours} is a vector of logicals,
                    unit-specific parameters \eqn{\phi_i} are included}
              \item{\code{linear}}{a \code{logical} (or a vector of logicals) indicating wether a linear 
                    trend \eqn{\beta} (or a linear trend \eqn{\beta_i} for each unit)
                    is included}
              \item{\code{nseason}}{Integer number of Fourier frequencies; if \code{nseason} is a vector
                    of integers, each unit \eqn{i} gets its own seasonal parameters }
              \item{\code{negbin}}{if \code{"single"} negative binomial rather than poisson is used,
                    if \code{"multiple"} unit-specific overdispersion parameters are used.}
              \item{\code{proportion}}{see details in \code{\link{meanResponse}} }
              }}        
        \item{thetastart}{vector with starting values for all parameters specified
            in the control object (for \code{optim}).  } 
        \item{verbose}{if \code{true} information about convergence is printed}

}
\value{Returns an object of class \code{ah} with elements 
   \item{coefficients}{estimated parameters}
   \item{se}{estimated standard errors}
   \item{cov}{covariance matrix}
   \item{loglikelihood}{loglikelihood}
   \item{convergence}{logical indicating whether \code{optim} converged or not}
   \item{fitted.values}{fitted mean values \eqn{\mu_{i,t}}{\mu_it} }
   \item{control}{specified control object}
   \item{disProgObj}{specified \code{disProg}-object}
}
\details{Note that for the time being this function is not a surveillance algorithm, but only a modelling approach as described in the Held et. al (2005) paper.}
\seealso{\code{\link{meanResponse}}}
\author{M. Paul, L. Held, M. Höhle}
\examples{

# univariate time series: salmonella agona cases 
data(salmonella.agona)
salmonella <- create.disProg(week=1:length(salmonella.agona$observed),
                             observed=salmonella.agona$observed, 
                             state=salmonella.agona$state)

model1 <- list(lambda=TRUE, linear=TRUE, 
               nseason=1, negbin="single")
					
algo.hhh(salmonella, control=model1)

# multivariate time series: 
# measles cases in Lower Saxony, Germany
data(measles.weser)

# same model as above
algo.hhh(measles.weser, control=model1)

# different starting values for 
# theta = (lambda, beta, gamma_1, gamma_2, psi)
startValues <- c(0.1, rep(0, 3), 1)
algo.hhh(measles.weser, control=model1, 
                 thetastart=startValues)

# include autoregressive parameter phi for adjacent "Kreise"
model2 <- list(lambda=TRUE, neighbours=TRUE, 
               linear=FALSE, nseason=1, 
               negbin="single")
algo.hhh(measles.weser, control=model2)
                  
## weekly counts of influenza and meningococcal infections in Germany, 2001-2006
data(influMen)

# specify model with two autoregressive parameters lambda_i, overdispersion
# parameters psi_i, an autoregressive parameter phi for meningococcal infections
# (i.e. nu_flu,t = lambda_flu * y_flu,t-1  
#  and  nu_men,t = lambda_men * y_men,t-1 + phi_men*y_flu,t-1 )
# and S=(3,1) Fourier frequencies
model <- list(lambda=c(TRUE,TRUE), neighbours=c(FALSE,TRUE),
              linear=FALSE,nseason=c(3,1),negbin="multiple")
              
# run algo.hhh
algo.hhh(influMen, control=model)

}
\keyword{ts}
\keyword{regression}
\source{Held, L., Höhle, M., Hofmann, M. (2005) A statistical framework for the 
    analysis of multivariate infectious disease surveillance counts. 
    \emph{Statistical Modelling}, \bold{5}, p. 187--199.}



% algo.hhh.grid
