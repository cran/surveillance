\name{anscombe.residuals}
\alias{anscombe.residuals}
\title{Compute Anscombe residuals}
\description{
  The residuals of \code{m} are transformed to form Anscombe residuals.
  which makes them approximately standard normal distributed.
}
\usage{
anscombe.residuals(m, phi)
}
\arguments{
  \item{m}{\code{m} is a glm object of the fit }
  \item{phi}{\code{phi} is the current estimated over-dispersion}
}
\value{Standardized Anscombe residuals of \code{m}}
\references{McCullagh & Nelder, Generalized Linear Models, 1989}
\keyword{}



%algo.farrington.assign.weights


