\name{algo.farrington.fitGLM}
\alias{algo.farrington.fitGLM}
\alias{algo.farrington.fitGLM.fast}
\alias{algo.farrington.fitGLM.populationOffset}
\title{Fit Poisson GLM of the Farrington procedure for a single time point}
\description{
  The function fits a Poisson regression model (GLM) with mean predictor
  \deqn{\log \mu_t = \alpha + \beta t}{
        log mu_t = alpha + beta * t}
  as specified by the Farrington procedure. If 
  requested, Anscombe residuals are computed based on an initial fit
  and a 2nd fit is made using weights, where base counts suspected to
  be caused by earlier outbreaks are downweighted.
}
\usage{
algo.farrington.fitGLM(response, wtime, timeTrend = TRUE,
                       reweight = TRUE, ...)
algo.farrington.fitGLM.fast(response, wtime, timeTrend = TRUE,
                            reweight = TRUE, ...)
algo.farrington.fitGLM.populationOffset(response, wtime, population,
                                        timeTrend=TRUE,reweight=TRUE, ...)
}
\arguments{
\item{response}{The vector of observed base counts}
\item{wtime}{Vector of week numbers corresponding to \code{response}}
\item{timeTrend}{Boolean whether to fit the \eqn{\beta t}{beta*t} or not}
\item{reweight}{Fit twice -- 2nd time with Anscombe residuals}
\item{population}{Population size. Possibly used as offset, i.e. in
    \code{algo.farrington.fitGLM.populationOffset} the value 
    \code{log(population)} is used as offset in the linear
    predictor of the GLM:   
    \deqn{\log \mu_t = \log(\texttt{population}) + \alpha + \beta t}{
      log mu_t = log(population) alpha + beta * t}
    This provides a way to adjust the Farrington procedure to the case
    of greatly varying populations. Note: This is an experimental implementation with methodology not covered by the original paper.
  }
\item{\dots}{Used to catch additional arguments, currently not used.}
}
\details{
  Compute weights from an initial fit and rescale using
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
\value{
  an object of class GLM with additional fields \code{wtime},
  \code{response} and \code{phi}. If the \code{glm} returns without
  convergence \code{NULL} is returned.
}
\seealso{\code{\link{anscombe.residuals}},\code{\link{algo.farrington}}}
\keyword{regression}
