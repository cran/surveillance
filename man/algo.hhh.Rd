% algo.hhh
\name{algo.hhh}
\alias{algo.hhh}
\alias{print.ah}
\alias{coef.ah}
\encoding{latin1}
\title{Model fit based on the Held, Hoehle, Hofman paper} 
\description{
        Fits a Poisson/negative binomial model with mean \eqn{\mu_{it}}{\mu_it} 
        (as described in Held/\enc{Höhle}{Hoehle}/Hofmann, 2005) to a multivariate 
        time series of counts.  
}
\usage{
algo.hhh(disProgObj, control=list(lambda=TRUE, neighbours=FALSE, 
   linear=FALSE, nseason = 0,
   negbin=c("none", "single", "multiple"), 
   proportion=c("none", "single", "multiple"),lag.range=NULL), 
   thetastart=NULL, verbose=TRUE)
}
\arguments{
  \item{disProgObj}{object of class \code{disProg}}
  \item{control}{control object:
        \describe{
        \item{\code{lambda}}{If \code{TRUE} an autoregressive 
              parameter \eqn{\lambda} is included, if \code{lambda} is a vector of logicals,
              unit-specific parameters \eqn{\lambda_i} are included. By default, observations \eqn{y_{t-lag}}{y_t-lag}
              at the previous time points, i.e. \eqn{lag=1}, are used for the autoregression. Other lags can be
              used by specifying \code{lambda} as a vector of integers, see Examples and Details.}
        \item{\code{neighbours}}{If \code{TRUE} an autoregressive parameter for 
              adjacent units \eqn{\phi} is included, if \code{neighbours} is a vector of logicals,
              unit-specific parameters \eqn{\phi_i} are included. By default, observations \eqn{y_{t-lag}}{y_t-lag}
              at the previous time points, i.e. \eqn{lag=1}, are used for the autoregression. Other lags can be
              used by specifying \code{neighbours} as a vector of integers.}
        \item{\code{linear}}{a \code{logical} (or a vector of logicals) indicating wether a linear 
              trend \eqn{\beta} (or a linear trend \eqn{\beta_i} for each unit)
              is included}
        \item{\code{nseason}}{Integer number of Fourier frequencies; if \code{nseason} is a vector
              of integers, each unit \eqn{i} gets its own seasonal parameters }
        \item{\code{negbin}}{if \code{"single"} negative binomial rather than poisson is used,
              if \code{"multiple"} unit-specific overdispersion parameters are used.}
        \item{\code{proportion}}{see Details}
        \item{\code{lag.range}}{determines which observations are used to fit the model }
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
    \item{lag}{which lag was used for the autoregressive parameters \eqn{lambda} and \eqn{phi} }
    \item{nObs}{number of observations used for fitting the model}
}
\details{
  This functions fits a model as specified in equations (1.2) and (1.1) in 
  Held et al. (2005) to univariate time series, and as specified in equations
  (3.3) and (3.2) (with extensions given in equations (2) and (4) in Paul et 
  al., 2008) to multivariate time series.

  For univariate time series, the mean structure of a Poisson or a negative
  binomial model is
   \deqn{\mu_t = \lambda y_{t-lag} + \nu_t }{\mu_t = \lambda y_t-lag + \nu_t }
   where  
   \deqn{\log( \nu_t) = \alpha + \beta t + \sum_{j=1}^{S}(\gamma_{2j-1} \sin(\omega_j t) + \gamma_{2j} \cos(\omega_j t) ) }{
         log(\nu_t) = \alpha + \beta t + 
                      \sum_(j=1)^S (\gamma_(2j-1) * sin(\omega_j * t) +
                                    \gamma_2j * cos(\omega_j * t) )  }
   and \eqn{\omega_j = 2\pi j/period }{\omega_j = 2 * \pi * j / period} are Fourier frequencies with
   known period, e.g. \code{period}=52 for weekly data.
   
   Per default, the number of cases at time point \eqn{t-1}, i.e. \eqn{lag=1}, enter
   as autoregressive covariates into the model. Other lags can also be considered.
     
   For multivariate time series the mean structure is
   \deqn{\mu_{it} = \lambda_i y_{i,t-lag} + \phi_i \sum_{j \sim i} w_{ji} y_{j,t-lag} + n_{it} \nu_{it}}{%
        \mu_it = \lambda_i * y_i,t-lag + 
                 \phi_i * \sum_(j ~ i) w_ji * y_j,t-lag + 
                 n_it * \nu_it }
   where 
   \deqn{\log(\nu_{it}) = \alpha_i + \beta_i t + \sum_{j=1}^{S_i} (\gamma_{i,2j-1} \sin(\omega_j t) + \gamma_{i,2j} \cos(\omega_j t) ) }{%
      log(\nu_it) = \alpha_i + \beta_i * t + 
                    \sum_(j=1)^S_i (\gamma_(i,2j-1) * sin(\omega_j * t) +
                                    \gamma_(i,2j) * cos(\omega_j * t) )  }
   and \eqn{n_{it}}{n_it} are standardized population counts. The weights \eqn{w_{ji}}{w_ji} are specified in the columns of
   the neighbourhood matrix \code{disProgObj$neighbourhood}.
   
   Alternatively, the mean can be specified as
   \deqn{\mu_{it} = \lambda_i \pi_i y_{i,t-1} + \sum_{j \sim i} \lambda_j (1-\pi_j)/ |k \sim j| y_{j,t-1} + n_{it} \nu_{it}}{%
       \mu_it = \lambda_i *\pi_i * y_i,t-1 + 
                \sum_(j ~ i) \lambda_j *(1-\pi_j)/|k ~ j| * y_j,t-1 + 
                n_it * \nu_it }
   if \code{proportion}="single" ("multiple") in the \code{control}
   argument. Note that this model specification is still experimental.   
}

\note{
For the time being this function is not a surveillance algorithm, 
but only a modelling approach as described in the papers by Held et. al (2005)  
and Paul et. al (2008).

}
\seealso{\code{\link{algo.hhh.grid}}, \code{\link{hhh4}}}
\author{M. Paul, L. Held, M. \enc{Höhle}{Hoehle}}
\examples{

# univariate time series: salmonella agona cases 
data(salmonella.agona)

model1 <- list(lambda=TRUE, linear=TRUE, 
               nseason=1, negbin="single")
               
algo.hhh(salmonella.agona, control=model1)

# multivariate time series: 
# measles cases in Lower Saxony, Germany
data(measles.weser)

# same model as above
algo.hhh(measles.weser, control=model1)

# include autoregressive parameter phi for adjacent "Kreise"
# specifiy start values for theta
model2 <- list(lambda = TRUE, neighbours = TRUE, 
               linear = FALSE, nseason = 1, 
               negbin = "single")
algo.hhh(measles.weser, control = model2, thetastart = rep(0, 20) )
                  
## weekly counts of influenza and meningococcal infections 
## in Germany, 2001-2006
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

# now meningococcal infections in the same week should enter as covariates
# (i.e. nu_flu,t = lambda_flu * y_flu,t-1  
#  and  nu_men,t = lambda_men * y_men,t-1 + phi_men*y_flu,t )
model2 <- list(lambda=c(1,1), neighbours=c(NA,0),
              linear=FALSE,nseason=c(3,1),negbin="multiple")
              
algo.hhh(influMen, control=model2)


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


