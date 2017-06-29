\name{estimateGLRNbHook}
\alias{estimateGLRNbHook}
\encoding{latin1}
\title{Hook function for in-control mean estimation}
\description{
 Estimation routine for the in-control mean of \code{\link{algo.glrpois}}.
 
 In \R < 2.14.0 and \pkg{surveillance} < 1.4 (i.e., without a package
 namespace) users could customize this function simply by
 defining a modified version in their workspace.
 This is no longer supported.
}
\usage{
estimateGLRNbHook()
}
\value{
  A list with elements
  \item{\code{mod}}{resulting model of a call of \code{glm.nb}}
  \item{\code{range}}{vector of length as \code{range} containing the predicted values}
}
\seealso{
\code{\link{algo.glrnb}}
}
\author{M. Hoehle}
\keyword{internal}
