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
  A list
  	\item{\code{mod}}{resulting model of a call of \code{glm.nb}}
    \item{\code{range}}{vector of length as \code{range} containing the predicted values} 

}
\seealso{
\code{\link{algo.glrnb}}
}
\author{M. Hoehle}
\examples{
\dontrun{
estimateGLRNbHook <- function() {
  #Fetch control object from parent
  control <- parent.frame()$control
  #The period
  p <- parent.frame()$disProgObj$freq
  #Current range to perform surveillance on
  range <- parent.frame()$range

  #Define training & test data set (the rest)
  train <- 1:(range[1]-1)
  test <- range
  
  #Perform an estimation based on all observations before timePoint
  #Event better - don't do this at all in the algorithm - force
  #user to do it himself - coz its a model selection problem
  data <- data.frame(y=parent.frame()$disProgObj$observed[t],t=train)
  #Build the model equation
  formula <- "y ~ 1 "
  if (control$mu0Model$trend) { formula <- paste(formula," + t",sep="") }
  for (s in 1:control$mu0Model$S) {
    formula <- paste(formula,"+cos(2*",s,"*pi/p*t)+ sin(2*",s,"*pi/p*t)",sep="")
  }
  #Fit the GLM
  m <- eval(substitute(glm.nb(form,data=data),
                      list(form=as.formula(formula))))

  #Predict mu_{0,t}
  return(as.numeric(predict(m,newdata=data.frame(t=range),type="response")))
}
}
}

\keyword{regression}
\keyword{internal}
