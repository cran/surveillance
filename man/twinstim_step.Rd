\name{twinstim_step}
\alias{stepComponent}
\alias{add1.twinstim}
\alias{drop1.twinstim}

\title{
  Stepwise Model Selection by AIC
}

\description{
\code{stepComponent} is a wrapper around \code{\link{step}} to select a
\code{"\link{twinstim}"} component's model based on an information
criterion in a stepwise algorithm.

There are also stand-alone single-step methods of \code{\link{add1}} and
\code{\link{drop1}}.
}

\usage{
stepComponent(object, component = c("endemic", "epidemic"),
              scope = list(upper = object$formula[[component]]),
              direction = "both", trace = 2, verbose = FALSE, ...)

\method{add1}{twinstim}(object, scope, component = c("endemic", "epidemic"), 
    trace = 2, ...)
\method{drop1}{twinstim}(object, scope, component = c("endemic", "epidemic"), 
     trace = 2, ...)
}

\arguments{
  \item{object}{an object of class \code{"twinstim"}.}
  \item{component}{one of \code{"endemic"} or \code{"epidemic"}
    (partially matched), determining the model component where the
    algorithm should proceed.}
  \item{scope,direction,trace}{see \code{\link{step}} and
    \code{\link{add1}}, respectively.}
  \item{verbose}{see \code{\link{twinstim}}.}
  \item{\dots}{further arguments passed to \code{\link{step}},
    \code{\link{add1.default}}, or \code{\link{drop1.default}},
    respectively.}
}

\value{
  See \code{\link{step}} and \code{\link{add1}}, respectively.
}

\author{
  (of this wrapper around \code{\link{step}}) Sebastian Meyer
}

\seealso{
\code{\link{step}}, \code{\link{add1}}, \code{\link{drop1}}
}

\examples{
data("imdepi", "imdepifit")

## simple baseline model
m0 <- update(imdepifit, epidemic=~1, siaf=NULL)

## AIC-based step-wise backward selection of the endemic component
m0_step <- stepComponent(m0, "endemic", scope=list(lower=~I(start/365-3.5)))
## nothing is dropped from the model

\dontshow{
m0_step$anova <- NULL
stopifnot(identical(m0, m0_step))
}
}

\keyword{models}
\keyword{methods}
