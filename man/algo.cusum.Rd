\name{algo.cusum}
\alias{algo.cusum}
\title{CUSUM method}
\description{ 
        Approximate one-side CUSUM method for a Poisson variate based on the 
        cumulative sum of the deviation between a reference value k and the 
        (standardized) observed values.
        An alarm is raised if the cumulative sum equals or exceeds a prespecified
        decision boundary h.
}
\usage{

  algo.cusum(disProgObj, control = list(range = range, k = 1.04, h = 2.26, 
             m = NULL, trans = "standard", alpha = NULL))
}
\arguments{
        \item{disProgObj}{object of class disProg (including the observed and the state chain).}
        \item{control}{control object: 
            \describe{
            \item{\code{range}}{determines the desired time points which should be evaluated}
            \item{\code{k}}{is the reference value}
            \item{\code{h}}{the decision boundary}
            \item{\code{m}}{expected number of cases -- if \code{NULL} then estimated from the data}
            \item{\code{trans}}{one of the following transformations (warning: anscombe and negbin transformations are experimental)
                \describe{
                \item{\code{rossi}}{compute standardized variables z3 as proposed by Rossi}
                \item{\code{standard}}{compute standardized variables z1 (based on asympotic normality)}
                \item{\code{anscombe}}{anscombe residuals -- experimental}
                \item{\code{anscombe2nd}}{ anscombe residuals as in Pierce and Schafer (1986) based on 2nd order approximation of E(X)  -- experimental}
                \item{\code{pearsonNegBin}}{compute Pearson residuals for NegBin -- experimental}
                \item{\code{anscombeNegBin}}{anscombe residuals for NegBin -- experimental}
                \item{\code{none}}{ no transformation}
                }
              }
              \item{\code{alpha}}{parameter of the negative binomial distribution, s.t. the variance is \eqn{m+alpha*m^2}}
        }
      }
}
\value{
\item{survRes}{
    \code{algo.cusum} gives a list of class \code{survRes} which includes the 
    vector of alarm values for every timepoint in \code{range} and the vector 
    of cumulative sums for every timepoint in \code{range} for the system 
    specified by \code{k} and \code{h}, the range and the input object of 
    class disProg.
  }
}
\details{The implementation is still experimental}
\seealso{}
\author{}
\examples{

    # Xi ~ Po(5), i=1,...,500
    disProgObj <- create.disProg(week=1:500, observed= rpois(500,lambda=5),
                                    state=rep(0,500))
    # there should be no alarms as mean doesn't change
    res <- algo.cusum(disProgObj, control = list(range = 100:500,trans="anscombe"))
    plot(res)
        
    # simulated data 
    disProgObj <- sim.pointSource(p = 1, r = 1, length = 250,
                              A = 0, alpha = log(5), beta = 0, phi = 10,
                              frequency = 10, state = NULL, K = 0)                           
    plot(disProgObj)
    
    # Test week 200 to 250 for outbreaks
    surv <- algo.cusum(disProgObj, control = list(range = 200:250))
    plot(surv)
}
\references{
  G. Rossi, L. Lampugnani and M. Marchi (1999), An approximate CUSUM procedure for surveillance of health events, Statistics in Medicine, 18, 2111--2122

  D. A. Pierce and D. W. Schafer (1986), Residuals in Generalized Linear Models, Journal of the American Statistical Association, 81, 977--986
}

\keyword{}
\source{}




