\name{hhh4_W_utils}
\alias{getNEweights}
\alias{coefW}
\title{
  Extract Neighbourhood Weights from a Fitted \code{hhh4} Model
}
\description{
  The \code{getNEweights} function extracts the (fitted) weight
  matrix/array from a \code{"hhh4"} object, after scaling and
  normalization.
  The \code{coefW} function extracts the coefficients of parametric
  neighbourhood weights from a \code{hhh4} fit (or directly from a
  corresponding coefficient vector), i.e., coefficients whose names
  begin with \dQuote{neweights}.
}
\usage{
getNEweights(object, pars = coefW(object),
             scale = ne$scale, normalize = ne$normalize)
coefW(object)
}
\arguments{
  \item{object}{an object of class \code{"hhh4"}.
    \code{coefW} also works with the coefficient vector.}
  \item{pars}{coefficients for parametric neighbourhood weights,
    such as for models using \code{\link{W_powerlaw}}.
    Defaults to the corresponding point estimates in \code{object}.}
  \item{scale,normalize}{parameters of the \code{ne} component of
    \code{\link{hhh4}}.}
}
\author{
  Sebastian Meyer
}
\keyword{utilities}
