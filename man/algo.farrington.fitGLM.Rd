\name{algo.farrington.fitGLM}
\alias{algo.farrington.fitGLM}
\alias{algo.farrington.fitGLM.fast}
\encoding{latin1}

\title{Fit the Poisson GLM of the Farrington procedure for a single
  time point}
\description{
  The function fits a Poisson regression model (GLM) with mean predictor
  \deqn{\log \mu_t = \alpha + \beta w_t}{
        log mu_t = alpha + beta * w}
  as specified by the Farrington procedure. That way we are able to
  predict the value \eqn{c_0}{c0}. If 
  requested Anscombe residuals are computed based on an initial fit
  and a 2nd fit is made using weights, where base counts suspected to
  be caused by earlier outbreaks are downweighted.
}
\usage{
  algo.farrington.fitGLM(response, wtime, timeTrend = TRUE, 
                         reweight = TRUE)
  algo.farrington.fitGLM.fast(response, wtime, timeTrend = TRUE, 
                         reweight = TRUE)
}
\arguments{
  \item{response}{The vector of observed base counts}
  \item{wtime}{Vector of week numbers corresponding to \code{response}}
  \item{timeTrend}{Boolean whether to fit the \eqn{\beta t}{beta*t} or not}
  \item{reweight}{Fit twice -- 2nd time with Anscombe residuals}
}
%
\details{Compute weights from an initial fit and rescale using
  Anscombe based residuals as described in the
  \code{\link{anscombe.residuals}} function.

  Note that \code{algo.farrington.fitGLM} uses the \code{glm} routine
  for fitting. A faster alternative is provided by
  \code{algo.farrington.fitGLM.fast} which uses the \code{glm.fit}
  function directly (thanks to Mikko Virtanen). This saves
  computational overhead and increases speed for 500 monitored time
  points by a factor of approximately two. However, some of the
  routine \code{glm} functions might not work on the output of this
  function. Which function is used for \code{algo.farrington} can be
  controlled by the \code{control$fitFun} argument.
  

}
%
\value{An object of class GLM with additional fields \code{wtime},
  \code{response} and \code{phi}. If the \code{glm} returns without
  convergence \code{NULL} is returned.}
%
\seealso{\code{\link{anscombe.residuals}},\code{\link{algo.farrington}}}
\keyword{regression}


