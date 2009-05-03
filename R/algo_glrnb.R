###################################################
### chunk number 1: 
###################################################

######################################################################
#
# Implementation of GLR -- documentation converted to Rd format.
#
# Author: Michael Hoehle (with contributions by Valentin Wimmer)
# Date:   8 Jan 2008
#
######################################################################

algo.glrnb <- function(disProgObj,
                       control = list(range=range,c.ARL=5,
                         mu0=NULL, alpha=0, Mtilde=1, M=-1, change="intercept",
                         theta=NULL,dir=c("inc","dec"),
                         ret=c("cases","value"))) {


  #Small helper function
  either <- function(cond, whenTrue, whenFalse) { if (cond) return(whenTrue) else return(whenFalse) }
  
  # Set the default values if not yet set
  if(is.null(control[["c.ARL",exact=TRUE]]))
    control$c.ARL <- 5
  if(is.null(control[["change",exact=TRUE]]))
    control$change <- "intercept"
  if(is.null(control[["Mtilde",exact=TRUE]]))
    control$Mtilde <- 1
  if(is.null(control[["M",exact=TRUE]]))
    control$M <- -1
  if(is.null(control[["dir",exact=TRUE]]))
    control$dir <- "inc"
  if(is.null(control[["ret",exact=TRUE]]))
  	control$ret <- "value"
  #if(is.null(control[["alpha",exact=TRUE]]))
  #    control$alpha <- 0

  #GLM (only filled if estimated)
  m <- NULL

  #Extract the important parts from the arguments
  observed <- disProgObj$observed
  t <- control$range
  control$mu0Model <- NULL
  range <- control$range
  control$dir <- match.arg(control$dir, c("inc","dec"))
  dir <- ifelse(control$dir=="inc",1,-1)
  control$ret <- match.arg(control$ret, c("value","cases"))
  ret <- pmatch(control$ret,c("value","cases"))
  mod <- list()
  
 

  # Estimate m (the expected number of cases), i.e. parameter lambda of a
  # poisson distribution based on time points 1:t-1
  if (is.null(control[["mu0",exact=TRUE]]) | is.list(control[["mu0",exact=TRUE]])) {
    #Initialize
    if (is.null(control[["mu0",exact=TRUE]])) control$mu0 <- list()
    if (is.null(control[["mu0",exact=TRUE]][["S"]])) control$mu0$S <- 1
    if (is.null(control[["mu0",exact=TRUE]][["trend"]])) control$mu0$trend <- FALSE
    if (is.null(control[["mu0",exact=TRUE]][["refit"]])) control$mu0$refit <- FALSE
    control$mu0Model <- control$mu0

    #Estimate using a hook function (lazy evaluation)
    control$mu0 <- estimateGLRNbHook()$pred
    
    mod[[1]] <- estimateGLRNbHook()$mod
    
    # if it is necessary to estimate alpha
    if(is.null(control[["alpha",exact=TRUE]])) control$alpha <- mod[[1]]$theta
  }
  
   #Postprocess
  if ((control$alpha>0) & (control$ret == "cases") & (is.null(control[["theta",exact=TRUE]]))) {
    stop("Return of cases is currently not implemented for the GLR detector based on the negative binomial distribution!")
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
    # cat("Doneidx === ",doneidx,"\n")
    # Call the C-interface -- this should depend on the type
    if (control$change == "intercept") {
      if (is.null(control[["theta",exact=TRUE]])) {
        if (control$alpha == 0) { #poisson

          if (control$M > 0 ){ # window limited
          
          	res <- .C("glr_cusum_window",as.integer(x),as.double(mu0),length(x),as.integer(control$M),as.integer(control$Mtilde),as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(dir),as.integer(ret),PACKAGE="surveillance")
        } 
        	else { # standard
        
        	res <- .C("glr_cusum",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(dir),as.integer(ret),PACKAGE="surveillance")

        	}
        } else { #negbin
          res <- .C("glr_nb_window",x=as.integer(x),mu0=as.double(mu0),alpha=as.double(control$alpha),lx=length(x),Mtilde=as.integer(control$Mtilde),M=as.integer(control$M),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),dir=as.integer(dir),PACKAGE="surveillance")
        }
      } else { ###################### !is.null(control$theta)
        if (control$alpha == 0) { #poisson

          res <- .C("lr_cusum",x=as.integer(x),mu0=as.double(mu0),lx=length(x),as.double(control$theta),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(ret),PACKAGE="surveillance")

        } else { #negbin
          warning("LR feature of the negative binomial distribution is currently experimental!")
          res <- .C("lr_cusum_nb",x=as.integer(x),mu0=as.double(mu0),alpha=as.double(control$alpha),lx=length(x),as.double(control$theta),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),cases=as.double(numeric(length(x))),as.integer(ret),PACKAGE="surveillance")

        }
      }
    } else { ################### Epidemic chart #######################
      if (control$change == "epi") {
        if (control$alpha == 0) { #pois
          res <- .C("glr_epi_window",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.integer(control$M),as.double(xm10),as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),PACKAGE="surveillance")
        } else {
          res <- .C("glr_nbgeneral_window",as.integer(x),as.double(mu0),alpha=as.double(control$alpha),lx=length(x),Mtilde=as.integer(control$Mtilde),M=as.integer(control$M),xm10=as.double(xm10),c.ARL=as.double(control$c.ARL),N=as.integer(0),val=as.double(numeric(length(x))),dir=as.integer(dir),PACKAGE="surveillance")
        }
      }
    }
    
     #In case an alarm found log this and reset the chart at res$N+1
    if (res$N < length(x)) {
      #Put appropriate value in upperbound
      upperbound[1:res$N + doneidx]  <- either(ret == 1, res$val[1:res$N] ,res$cases[1:res$N])
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
        mu0 <- estimateGLRNbHook()$pred
        mod[[noofalarms+2]] <-  estimateGLRNbHook()$mod 
        control$mu0[(doneidx + res$N + 1):length(control$mu0)] <- mu0
      }

      noofalarms <- noofalarms + 1

    }
    doneidx <- doneidx + res$N
  }



  # fix of the problem that no upperbound-statistic is returned after 
  #last alarm
  upperbound[(doneidx-res$N+1):nrow(upperbound)] <- either(ret == 1, res$val, res$cases)
  
  #fix of the problem that no upperbound-statistic ss returned 
  #in case of no alarm
  if (noofalarms == 0) {
    upperbound <- either(ret==1, res$val, res$cases)
  }

  # ensure upper bound is positive and not NaN
  upperbound[is.na(upperbound)] <- 0
  upperbound[upperbound < 0] <- 0
  
 
  # Add name and data name to control object
  
  algoName <- either(control$alpha == 0, "glrpois:", "glrnb:")
  control$name <- paste(algoName, control$change)
  control$data <- paste(deparse(substitute(disProgObj)))
  control$m    <- m
  control$mu0Model$fitted <- mod

  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, 
                 disProgObj=disProgObj,control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}


###################################################
### chunk number 2: 
###################################################
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
  data <- data.frame(y=parent.frame()$disProgObj$observed[train],t=train)
  #Build the model equation
  formula <- "y ~ 1 "
  if (control$mu0Model$trend) { formula <- paste(formula," + t",sep="") }
  for (s in 1:control$mu0Model$S) {
    formula <- paste(formula,"+cos(2*",s,"*pi/p*t)+ sin(2*",s,"*pi/p*t)",sep="")
  }
  #Fit the GLM
  m <- eval(substitute(glm.nb(form,data=data),list(form=as.formula(formula))))

  #Predict mu_{0,t}
  return(list(mod=m,pred=as.numeric(predict(m,newdata=data.frame(t=range),type="response"))))
}


