\name{anscombe.residuals}
\alias{anscombe.residuals}

\title{Compute Anscombe Residuals}
\description{
  Compute Anscombe residuals from a fitted \code{\link{glm}},
  which makes them approximately standard normal distributed.
}
\usage{
anscombe.residuals(m, phi)
}
\arguments{
  \item{m}{a fitted \code{"glm"}}
  \item{phi}{the current estimated overdispersion}
}
\value{The standardized Anscombe residuals of \code{m}}
\references{McCullagh & Nelder, Generalized Linear Models, 1989}
\keyword{regression}
