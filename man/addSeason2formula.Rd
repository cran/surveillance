\name{addSeason2formula}
\alias{addSeason2formula}
\title{
  Add Harmonics to an Existing Formula
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
    \sin( s \cdot 2\pi \cdot \code{timevar}/\code{period} ),\;
    \cos( s \cdot 2\pi \cdot \code{timevar}/\code{period} ),
  }{
    sin( s * 2 * pi * \code{timevar}/\code{period} ),
    cos( s * 2 * pi * \code{timevar}/\code{period} ),
  }
  for \eqn{s = 1,\dots,\code{S}} to an existing formula \code{f}.

  Note the following equivalence when interpreting the coefficients of
  the seasonal terms:
  \deqn{
    \gamma \sin(\omega t) + \delta \cos(\omega t) = 
    A \sin(\omega t + \epsilon)
  }{
    gamma sin(omega * t) + delta cos(omega * t) = 
    A sin(omega * t + epsilon)
  }
  with amplitude \eqn{A = \sqrt{\gamma^2 + \delta^2}}{A = sqrt(gamma^2 + delta^2)}
  and phase shift \eqn{\epsilon = \arctan(\delta / \gamma)}. 
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
