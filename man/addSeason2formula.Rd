\name{addSeason2formula}
\alias{addSeason2formula}
\title{
  Function that adds a sine-/cosine formula to an existing formula.
}
\description{
  This function helps to construct a \code{\link{formula}} object that
  can be used in a call to \code{\link{hhh4}} to model 
  seasonal variation via a sum of sine and cosine terms. 
}
\usage{
addSeason2formula(f = ~1, S = 1, period = 52, timevar = "t")
}
\arguments{
  \item{f}{
  formula that the seasonal terms should be added to, 
  defaults to an intercept \code{~1}.
}
  \item{S}{
  number of sine and cosine terms. If \code{S} is a vector, 
  unit-specific seasonal terms are created.
  
}
  \item{period}{
  period of the season, defaults to 52 for weekly data.
}
  \item{timevar}{
    the time variable in the model. Defaults to \code{"t"}.
  }
}
\details{
  The function adds the seasonal terms
  \deqn{
        \sum_{s=1}^\code{S} \gamma_s \sin(\frac{2\pi s}{\code{period}} t)
                    +\delta_s \cos(\frac{2\pi s}{\code{period}} t),
        }{
        sum_s gamma_s * sin(2*pi*s/period * t) 
              + delta_s * cos(2*pi*s/period * t),
        }
  where \eqn{\gamma_s}{gamma_s} and \eqn{\delta_s}{delta_s} are the unknown 
  parameters and \eqn{t}, \eqn{t = 1, 2, \ldots} denotes the time
  variable \code{timevar}, to an existing formula \code{f}.

  Note that the seasonal terms can also be expressed as
  \deqn{\gamma_{s} \sin(\frac{2\pi s}{\code{period}} t) + \delta_{s} \cos(\frac{2\pi s}{\code{period}} t) = 
        A_s \sin(\frac{2\pi s}{\code{period}} t + \epsilon_s)}{%
      \gamma_s sin(2*pi*s/period * t) + \delta_s cos2*pi*s/period * t) = 
        A_s sin(2*pi*s/period * t + \epsilon_s)}
  with amplitude \eqn{A_s=\sqrt{\gamma_s^2 +\delta_s^2}}{A_s=sqrt{\gamma_s^2 +\delta_s^2}} 
  and phase shift \eqn{\tan(\epsilon_s) = \delta_s / \gamma_s}. 
  The amplitude and phase shift can be obtained from a fitted
  \code{\link{hhh4}} model via \code{coef(..., amplitudeShift = TRUE)},
  see \code{\link{coef.hhh4}}.
}
\value{
  Returns a \code{\link{formula}} with the seasonal terms added and 
  its environment set to \code{\link{.GlobalEnv}}.
  Note that to use the resulting formula in \code{\link{hhh4}}, 
  a time variable named as specified by the argument \code{timevar} must
  be available.
}

\author{
M. Paul, with contributions by S. Meyer
}

\seealso{
\code{\link{hhh4}}, \code{\link{fe}}, \code{\link{ri}}
}

\examples{
# add 2 sine/cosine terms to a model with intercept and linear trend
addSeason2formula(f = ~ 1 + t, S = 2)

# the same for monthly data
addSeason2formula(f = ~ 1 + t, S = 2, period = 12)

# different number of seasons for a bivariate time series
addSeason2formula(f = ~ 1, S = c(3, 1), period = 52)
}
