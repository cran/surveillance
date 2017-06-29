\name{sim.seasonalNoise}
\alias{sim.seasonalNoise}
\encoding{latin1}

\title{Generation of Background Noise for Simulated Timeseries}
\description{Generation of a cyclic model of a Poisson distribution
as background data for a simulated timevector.

The mean of the Poisson distribution is modelled as:
\deqn{\mu = \exp(A \sin( frequency \cdot \omega \cdot (t + \phi)) + \alpha  + \beta * t +  K * state)}{%
mu = exp(A * sin( frequency * omega * (t + phi)) + alpha  + beta * t +  K * state)}
}
\usage{
sim.seasonalNoise(A = 1, alpha = 1, beta = 0, phi = 0,
                  length, frequency = 1, state = NULL, K = 0)
}
\arguments{
  \item{A}{amplitude (range of sinus), default = 1.}
  \item{alpha}{parameter to move along the y-axis (negative values not allowed)
    with alpha > = A, default = 1.}
  \item{beta}{regression coefficient, default = 0.}
  \item{phi}{factor to create seasonal moves
    (moves the curve along the x-axis), default = 0.}
  \item{length}{number of weeks to model.}
  \item{frequency}{factor to determine the oscillation-frequency, default = 1.}
  \item{state}{if a state chain is entered the outbreaks will be additional
    weighted by K.}
  \item{K}{additional weigth for an outbreak which influences the distribution
    parameter mu, default = 0.}
}
\value{
  an object of class \code{seasonNoise} which includes the modelled
  timevector, the parameter \code{mu} and all input parameters.
}
\seealso{\code{\link{sim.pointSource}}}
\author{M. \enc{Höhle}{Hoehle}, A. Riebler, C. Lang}
\examples{
season <- sim.seasonalNoise(length = 300)
plot(season$seasonalBackground,type = "l")

# use a negative timetrend beta
season <- sim.seasonalNoise(beta = -0.003, length = 300)
plot(season$seasonalBackground,type = "l")
}
\keyword{datagen}
