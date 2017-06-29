\name{xtable.algoQV}
\alias{xtable.algoQV}
\title{Xtable quality value object}
\description{xtable a single quality value object in a nicely formatted way}
\usage{
\method{xtable}{algoQV}(x,caption = NULL, label = NULL, 
      align = NULL, digits = NULL, display = NULL, ...)
}
\arguments{
\item{x}{Quality Values object generated with \code{quality}}
\item{caption}{See \code{\link[xtable]{xtable}}}
\item{label}{See \code{\link[xtable]{xtable}}}
\item{align}{See \code{\link[xtable]{xtable}}}
\item{digits}{See \code{\link[xtable]{xtable}}}
\item{display}{See \code{\link[xtable]{xtable}}}
\item{...}{Further arguments (see \code{\link[xtable]{xtable})}}
}
\keyword{print}
\seealso{ \code{\link[xtable]{xtable}}}
\examples{
# Create a test object
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 200, A = 1,
                              alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.7)

# Let this object be tested from rki1
survResObj <- algo.rki1(disProgObj, control = list(range = 50:200))

# Compute the quality values in a nice formatted way
xtable(algo.quality(survResObj))
}
