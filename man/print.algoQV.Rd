\name{print.algoQV}
\alias{print.algoQV}
\title{Print Quality Value Object}
\description{Print a single quality value object in a nicely formatted way}
\usage{
  \method{print}{algoQV}(x,...)
}
\arguments{
\item{x}{Quality Values object generated with \code{quality}}
\item{...}{Further arguments (not really used)}
}
\examples{
# Create a test object
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 200, A = 1,
                              alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.7)

# Let this object be tested from rki1
survResObj <- algo.rki1(disProgObj, control = list(range = 50:200))

# Compute the quality values in a nice formatted way
algo.quality(survResObj) 
}
\keyword{print}
