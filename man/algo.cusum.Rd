\name{algo.cusum}
\alias{algo.cusum}
\title{CUSUM method}
\encoding{latin1}
\description{ 
        Approximate one-side CUSUM method for a Poisson variate based on the 
        cumulative sum of the deviation between a reference value k and the 
        transformed observed values.
        An alarm is raised if the cumulative sum equals or exceeds a prespecified
        decision boundary h. The function can handle time varying expectations.
}
\usage{
algo.cusum(disProgObj, control = list(range = range, k = 1.04, h = 2.26, 
           m = NULL, trans = "standard", alpha = NULL))
}
\arguments{
        \item{disProgObj}{object of class disProg (including the observed and the state chain)}
        \item{control}{control object: 
            \describe{
            \item{\code{range}}{determines the desired time points which should be evaluated}
            \item{\code{k}}{is the reference value}
            \item{\code{h}}{the decision boundary}
            \item{\code{m}}{how to determine the expected number of cases -- 
                the following arguments are possible
                \describe{
                \item{\code{numeric}}{a vector of values having the
                    same length as \code{range}. If a single numeric
                    value is specified then this value is replicated
                    \code{length(range)} times.}
                \item{\code{NULL}}{A single value is estimated by
                    taking the mean of all observations previous to
                    the first \code{range} value.}
                \item{\code{"glm"}}{ A GLM of the form \deqn{\log(m_t)
                      = \alpha + \beta t + \sum_{s=1}^S (\gamma_s
                      \sin(\omega_s t) + \delta_s \cos(\omega_s t)),}
                    where \eqn{\omega_s = \frac{2\pi}{52}s}{\omega_s = 2\pi/52 s} are the
                    Fourier frequencies is fitted. Then this model is
                    used to predict the \code{range} values.}
                }}
            \item{\code{trans}}{one of the following transformations (warning: Anscombe and NegBin transformations are experimental)
                \describe{
                \item{\code{rossi}}{standardized variables z3 as proposed by Rossi}
                \item{\code{standard}}{standardized variables z1 (based
             on asymptotic normality) - This is the default.}
                \item{\code{anscombe}}{anscombe residuals -- experimental}
                \item{\code{anscombe2nd}}{ anscombe residuals as in Pierce and Schafer (1986) based on 2nd order approximation of E(X)  -- experimental}
                \item{\code{pearsonNegBin}}{compute Pearson residuals for NegBin -- experimental}
                \item{\code{anscombeNegBin}}{anscombe residuals for NegBin -- experimental}
                \item{\code{none}}{ no transformation}
                }
              }
              \item{\code{alpha}}{parameter of the negative binomial distribution, s.t. the variance is \eqn{m+\alpha *m^2} }
        }
      }
}
\value{
    \code{algo.cusum} gives a list of class \code{"survRes"} which includes the 
    vector of alarm values for every timepoint in \code{range} and the vector 
    of cumulative sums for every timepoint in \code{range} for the system 
    specified by \code{k} and \code{h}, the range and the input object of 
    class \code{"disProg"}.

    The \code{upperbound} entry shows for each time instance the number of diseased individuals
    it would have taken the cusum to signal. Once the CUSUM signals no resetting is applied, i.e.
    signals occurs until the CUSUM statistic again returns below the threshold.

    In case \code{control$m="glm"} was used, the returned
    \code{control$m.glm} entry contains the fitted \code{"glm"} object.
}
\note{This implementation is experimental, but will not be developed further.}
\author{M. Paul and M. \enc{Höhle}{Hoehle}}
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

\keyword{classif}
