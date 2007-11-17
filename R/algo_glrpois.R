###################################################
### chunk number 1: 
###################################################

######################################################################
# 
# Implementation of GLR -- documentation converted to Rd format.
#
# Author: Michael Hoehle
# Date:   27 Nov 2006
#
######################################################################

algo.glrpois <- function(disProgObj, 
                         control = list(range=range,c.ARL=5, 
                           mu0=NULL, Mtilde=1, M=-1, change="intercept",theta=NULL,dir="inc")){
  
  # Set the default values if not yet set
  if(is.null(control$c.ARL))
    control$c.ARL <- 5
  if(is.null(control$change))
    control$change <- "intercept" 
  if(is.null(control$Mtilde))
    control$Mtilde <- 1
  if(is.null(control$M))
    control$M <- -1
  if(is.null(control$dir))
    control$dir <- "inc"

  #GLM (only filled if estimated)
  m <- NULL

  #Extract the important parts from the arguments
  observed <- disProgObj$observed
  t <- control$range
  control$mu0Model <- NULL
  range <- control$range
  dir <- ifelse(control$dir=="inc",1,-1)

  # Estimate m (the expected number of cases), i.e. parameter lambda of a
  # poisson distribution based on time points 1:t-1
  if (is.null(control$mu0) | is.list(control$mu0)) {
    #Initialize
    if (is.null(control$mu0)) control$mu0 <- list()
    if (is.null(control$mu0$S)) control$mu0$S <- 1
    if (is.null(control$mu0$trend)) control$mu0$trend <- FALSE
    if (is.null(control$mu0$refit)) control$m0$refit <- FALSE
    control$mu0Model <- control$mu0

    #Estimate using a hook function (lazy evaluation)
    control$mu0 <- estimateGLRPoisHook()
  } 
  
  #The counts
  x <- observed[control$range]
  mu0 <- control$mu0

  #Reserve space for the results
  # start with cusum[timePoint -1] = 0, i.e. set cusum[1] = 0
  alarm <- matrix(data = 0, nrow = length(t), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(t), ncol = 1)

  #Setup counters for the progress
  doneidx <- 0
  N <- 1
  xm10 <- 0
  noofalarms <- 0
  noOfTimePoints <- length(t)
  #Loop as long as we are not through the sequence
  while (doneidx < noOfTimePoints) {
    #cat("Doneidx === ",doneidx,"\n")
    #Call the C-interface -- this should depend on the type
    if (control$change == "intercept") {
      if (is.null(control$theta)) {
        res <- .C("glr_cusum",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.double(control$c.ARL),N=as.integer(0),val=as.double(x),as.integer(dir),PACKAGE="surveillance")
      } else {
        res <- .C("lr_cusum",as.integer(x),as.double(mu0),length(x),as.double(control$theta),as.double(control$c.ARL),N=as.integer(0),val=as.double(x),PACKAGE="surveillance")
      }
    } else {
      ########################## Epidemic chart #######################
      if (control$change == "epi") {
        res <- .C("glr_epi_window",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.integer(control$M),as.double(xm10),as.double(control$c.ARL),N=as.integer(0),val=as.double(x),PACKAGE="surveillance")
      }
    }

    #In case an alarm found log this and reset the chart at res$N+1
    if (res$N < length(x)) {
      upperbound[1:res$N + doneidx] <- res$val[1:res$N]
      alarm[res$N + doneidx] <- TRUE

      #Chop & get ready for next round
      xm10 <- x[res$N] #put start value x_0 to last value
      x <- x[-(1:res$N)] ; t <- t[-(1:res$N)] 
      #If no refitting is to be done things are easy
      if (!is.list(control$mu0Model) || (control$mu0Model$refit == FALSE)) { 
        mu0 <- mu0[-(1:res$N)]
      } else {
        #Update the range (how to change back??)
        range <- range[-(1:res$N)]
        mu0 <- estimateGLRPoisHook()
        control$mu0[(doneidx + res$N + 1):length(control$mu0)] <- mu0
      }
      
			 noofalarms <- noofalarms + 1
      
    }
    doneidx <- doneidx + res$N
  }

	# fix of the problem that no upperbound-statistic was returned in case of no alarm
	if (noofalarms==0) upperbound <- res$val 

  # ensure upper bound is positive and not NaN
  upperbound[is.na(upperbound)] <- 0
  upperbound[upperbound < 0] <- 0
    
  
  #Add name and data name to control object.
  control$name <- paste("glrpois:", control$change)
  control$data <- paste(deparse(substitute(disProgObj)))
  control$m    <- m
  
  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj,control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}



###################################################
### chunk number 2: 
###################################################
estimateGLRPoisHook <- function() {
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
  m <- eval(substitute(glm(form,family=poisson(),data=data),list(form=as.formula(formula))))

  #Predict mu_{0,t}
  return(as.numeric(predict(m,newdata=data.frame(t=range),type="response")))
}


