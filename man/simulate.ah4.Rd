\encoding{latin1}
\name{simulate.ah4}
\alias{simulate.ah4}

\title{Simulates data based on the model proposed by Paul and Held (2011)
}

\description{
   Simulates a multivariate time series of counts based on the
   Poisson/Negative Binomial model as described in Paul and Held (2011).
}

\usage{
\method{simulate}{ah4}(object, nsim = 1, seed = NULL, y.start = NULL, coefs = NULL, ...)
}

\arguments{
  \item{object}{
  an object of class \code{"ah4"}.
}
  \item{nsim}{
  number of time series to simulate. Defaults to \code{1}.
}
  \item{seed}{
   an integer that will be used in the call to \code{\link{set.seed}} before
   simulating the time series.
}
  \item{y.start}{
   vector with starting counts for the epidemic components. 
   If \code{NULL}, the observed means in the respective units of the
   data in \code{object} are used.
}
  \item{coefs}{
   if not \code{NULL}, these values (in the same order as \code{coef(object)})
   are used to simulate from the model specified in \code{object}.
}
  \item{\dots}{unused (argument of the generic).
}
}

\details{
  Simulates data from a Poisson or a Negative Binomial model
  with mean
  \deqn{\mu_{it} = \lambda_{it} y_{i,t-1} + 
                   \phi_{it} \sum_{j \neq i} w_{ji} y_{j,t-1} + 
                   \nu_{it}}{%
        \mu_it = \lambda_it y_i,t-1 + 
                  \phi_it \sum_j w_ji y_j,t-1 + 
                  \nu_it}
  where
  \eqn{\lambda_{it}>0}, \eqn{\phi_{it}>0}, and \eqn{\nu_{it}>0} are 
  parameters which are modelled parametrically.
  The function uses the model and parameter estimates of the fitted
  \code{object} to simulate the time series.

  With the argument \code{coefs} it is possible to simulate from  
  the model as specified in \code{object}, but with different 
  parameter values.
}

\value{
An object of class \code{"sts"} in the case of \code{nsim = 1}, and
a list of \code{"sts"} objects otherwise.
}

\references{
  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136
}
\author{
M. Paul
}

\seealso{
\code{\link{hhh4}}, \code{\link{simHHH}}
}

\examples{
data(influMen)
# convert to sts class and extract meningococcal disease time series
meningo <- disProg2sts(influMen)[,2]

# fit model
fit <- hhh4(meningo, control = list(ar = list(f = ~ 1),
            end = list(f = addSeason2formula(S = 1, period = 52)),
            family = "NegBin1"))
plot(fit)

# simulate from model
set.seed(1234)
simData <- simulate(fit)

# plot simulated data
plot(simData, main = "simulated data", legend.opts = NULL, xaxis.years = FALSE)

# consider a Poisson instead of a NegBin model
coefs <- coef(fit)
coefs["overdisp"] <- 0

simData2 <- simulate(fit, coefs = coefs)
plot(simData2, main = "simulated data: Poisson model", 
     legend.opts = NULL, xaxis.years = FALSE)

# consider a model with higher autoregressive parameter
coefs <- coef(fit)
coefs[1] <- log(0.5)

simData3 <- simulate(fit, coefs = coefs)
plot(simData3, main = "simulated data: lambda = 0.5", 
     legend.opts = NULL, xaxis.years = FALSE)

}

\keyword{datagen}
