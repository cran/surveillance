\name{hhh4_predict}
\alias{predict.hhh4}
\title{Predictions from a \code{hhh4} Model}
\description{
  Get fitted (component) means from a \code{\link{hhh4}} model.
}

\usage{
\method{predict}{hhh4}(object, newSubset=object$control$subset,
        type="response", \dots)
}

\arguments{
\item{object}{fitted \code{\link{hhh4}} model (class \code{"hhh4"}).}
\item{newSubset}{subset of time points for which to return the
  predictions. Defaults to the subset used for fitting the model, and
  must be a subset of \code{1:nrow(object$stsObj)}.}
\item{type}{the type of prediction required. The default
  (\code{"response"} or, equivalently, \code{"mean"}) is on the 
  scale of the response variable (mean = endemic plus epidemic components).
  The alternatives are: \code{"endemic"}, \code{"epidemic"},
  \code{"epi.own"} (i.e. the autoregressive part), and
  \code{"epi.neighbours"} (i.e. the spatio-temporal part).}
\item{\dots}{unused (argument of the generic).}
}

\value{
  matrix of fitted means for each time point (of \code{newSubset}) and region.
}

\author{Michaela Paul and Sebastian Meyer}

\keyword{methods}
\keyword{models}
