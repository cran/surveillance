##############################################################################
# This function is a wrapper for univariate surveillance algorithms
# using the old disProg and survRes object
#
# An sts object is given and a pre specified algorithms is ran
# by successively creating a disProg object for each region,
# running the algo and then assign the slots of the resulting survRes
# object to an sts object.
###################################################################################

###Apply other algorithms by wrapping up a suitable package.

#Wrapper function to call algo.farrington for each time series in an sts object
wrap.algo <- function(sts, algo, control,
                      control.hook=function(k) return(control),...) {
  #Number of time series
  nAreas <- ncol(sts@observed)
  nAlarm <- length(control$range)
  #Adjust alarm matrix so we have the same number of values as control$range
  sts@alarm <- matrix(NA,ncol=nAreas,nrow=nAlarm)
  sts@upperbound <- matrix(NA,ncol=nAreas,nrow=nAlarm)

  #Loop
  for (k in 1:nAreas) {
    cat("Running ",algo," on area ",k," out of ",nAreas,"\n")
    
    ##Create an old S4 disProg object
    disProg.k <- create.disProg(sts@week, sts@observed[,k], sts@state[,k], freq=sts@freq, start=sts@start)

    #Use the univariate algorithm (possibly preprocess control object)
    kcontrol <- control.hook(k)
    survRes.k <- do.call(algo,args = list(disProg.k, control=kcontrol))
    
    #Transfer results to the S4 object
    if (!is.null(survRes.k)) {
      sts@alarm[,k] <- survRes.k$alarm
      sts@upperbound[,k] <- survRes.k$upperbound

      #Control object needs only to be set once
      sts@control <- survRes.k$control
    }
  }

  #Throw away the un-needed observations
  sts@observed <- sts@observed[control$range,,drop=FALSE]
  sts@state <- sts@state[control$range,,drop=FALSE]
  sts@populationFrac <- sts@populationFrac[control$range,,drop=FALSE]
  
  #Fix the corresponding start entry
  start <- sts@start
  new.sampleNo <- start[2] + min(control$range) - 1
  start.year <- start[1] + (new.sampleNo - 1) %/% sts@freq 
  start.sampleNo <- (new.sampleNo - 1) %% sts@freq + 1
  sts@start <- c(start.year,start.sampleNo)

  #Ensure dimnames
  sts <- fix.dimnames(sts)
  
  return(sts)
}

#Farrington wrapper
farrington <- function(sts, control=list(range=NULL, b=3, w=3, reweight=TRUE, verbose=FALSE,alpha=0.01)) {
  wrap.algo(sts,algo="algo.farrington",control=control)
}

#CDC wrapper
cdc <- function(sts, control= list(range = range,alpha=0.025)) {
  wrap.algo(sts,algo="algo.cdc",control=control)
}

#Bayes wrapper (this can be implemented more efficiently)
bayes <- function(sts, control = list(range = range, b = 0, w = 6, actY = TRUE,alpha=0.05)) {
  wrap.algo(sts,algo="algo.bayes",control=control)
}

#RKI wrapper
rki <- function(sts, control = list(range = range, b = 2, w = 4, actY = FALSE)) {
  wrap.algo(sts,algo="algo.rki",control=control)
}

#Cusum wrapper
cusum <- function(sts,  control = list(range=range, k=1.04, h=2.26, m=NULL, trans="standard",alpha=NULL)) {
  wrap.algo(sts,algo="algo.cusum",control=control)
}

#GLRpois wrapper
glrpois <- function(sts, control = list(range=range,c.ARL=5, S=1,
                           beta=NULL, Mtilde=1, M=-1, change="intercept",theta=NULL)) {
  wrap.algo(sts,algo="algo.glrpois",control=control)
}

#### this code definitely needs some more documentation -- wrap.algo atm is
# 100% without docu
#Rogerson wrapper
# theta0t now has to be a matrix
#library(surveillance)
#data("ha")
#rogerson(disProg2sts(ha),control=list(range=200:290,ARL0=100,s=1,theta0t=matrix(1,nrow=91,ncol=12)))

rogerson <- function(sts, control = list(range=range, theta0t=NULL,
                            ARL0=NULL, s=NULL, hValues=NULL,
                            distribution=c("poisson","binomial"),
                            nt=NULL, FIR=FALSE,limit=NULL, digits=1)) {
  #Hook function to find right theta0t vector
  control.hook = function(k) {
    control$hValues <- hValues(theta0 = control$theta0t[,k], ARL0=control$ARL0, control$s , distr = control$distribution)$hValues
    control$theta0t <- control$theta0t[,k]
#    print(control)
    return(control)
  }
  #WrapIt
  wrap.algo(sts,algo="algo.rogerson",control=control,control.hook=control.hook)
}


