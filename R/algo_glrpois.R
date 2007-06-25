###################################################
### chunk number 1: 
###################################################

######################################################################
# 
# Implementation of GLR -- documentation converted to Rd format.
#
# Author: Michael Höhle
# Date:   27 Nov 2006
#
######################################################################

algo.glrpois <- function(disProgObj, 
                         control = list(range=range,c.ARL=5, S=1,
                           mu0=NULL, Mtilde=1, M=-1, change="intercept",theta=NULL)){
  
  # Set the default values if not yet set
  if(is.null(control$c.ARL))
    control$c.ARL <- 5
  if(is.null(control$change))
    control$change <- "intercept" #nothing else handles atm.
  if(is.null(control$Mtilde))
    control$Mtilde <- 1
  if(is.null(control$M))
    control$M <- -1
  if(is.null(control$S)) {
    control$S <- 1
  } else {
    if(control$S>1) cat("Error: Only S=1 is handled atm.\n")
  }

  #Extract the important parts from the arguments
  observed <- disProgObj$observed
  timePoint <- control$range[1] #optimize on this
  t <- control$range
  #The period
  p <- disProgObj$freq

  # Estimate m (the expected number of cases), i.e. parameter lambda of a
  # poisson distribution based on time points 1:t-1
  if(is.null(control$mu)) {
    #Perform an estimation based on all observations before timePoint
    #NOT DONE YET -- if S>1
    t <- 1:(timePoint-1)
    data <- data.frame(x=disProgObj$observed[t],t=t)
    #different model!
    m <- glm(x ~ 1 + cos(2*pi/p*t)+ sin(2*pi/p*t),family=poisson(),data=data)
    beta <- coef(m)

    #Compute mu_{0,t}
    omega <- 2*pi/disProgObj$freq
    base <- rep(beta[1],length(t))
    base <- for (s in 1:control$S) {
      base <- base + beta[2*s] * cos(omega*t*s) + beta[2*s+1]*sin(omega*t*s)
    }
    mu0 <- exp(base)
  }
  
  #The counts
  x <- observed[control$range]

  #Reserve space for the results
  # start with cusum[timePoint -1] = 0, i.e. set cusum[1] = 0
  alarm <- matrix(data = 0, nrow = length(t), ncol = 1)

  #Call the C-interface -- this should depend on the type
  #Force computations until the end
  if (control$change == "intercept") {
    if (is.null(control$theta)) {
      res <- .C("glr_cusum",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.double(1e5),N=as.integer(0),val=as.double(x),PACKAGE="surveillance")
    } else {
      res <- .C("lr_cusum",as.integer(x),as.double(mu0),length(x),as.double(control$theta),as.double(1e5),N=as.integer(0),val=as.double(x),PACKAGE="surveillance")
    }
  } else {
    ########################## Epidemic chart #######################
    if (control$change == "epi") {
      res <- .C("glr_epi_window",as.integer(x),as.double(mu0),length(x),as.integer(control$Mtilde),as.integer(control$M),as.double(1e5),N=as.integer(0),val=as.double(x),PACKAGE="surveillance")
    }
  }
  
  #Retrieve results -- this should also include some val's!
  #print(res$N)
  #N <- res$N + min(t) -1 ; alarm[res$N] <- 1
  alarm <- res$val > control$c.ARL
  

  #upperbound <- rep(0,length(time))
  upperbound <- res$val
  # ensure upper bound is positive and not NaN
  upperbound[is.na(upperbound)] <- 0
  upperbound[upperbound < 0] <- 0
    
  
  #Add name and data name to control object.
  control$name <- paste("seaspois:", control$change)
  control$data <- paste(deparse(substitute(disProgObj)))
  control$m    <- m
  
  # return alarm and upperbound vectors
  result <- list(alarm = alarm, upperbound = upperbound, disProgObj=disProgObj,control=control)

  class(result) = "survRes" # for surveillance system result
  return(result)
}







