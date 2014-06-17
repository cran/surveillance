######################################################################
# An implementation of the Bayesian Outbreak Detection Algorithm (BODA)
# described in Manitz and H{\"o}hle (2013), Biometrical Journal.
#
# Note: The algorithm requires the non CRAN package INLA to run.
# You can install this package by writing
# R> source("http://www.math.ntnu.no/inla/givemeINLA.R")
# See http://www.r-inla.org/download for details.
#
#
# Author:
# The initial code was written by J. Manitz, which was then later
# adapted and modified for integration into the package by M. Hoehle.
#
# Date:
#  Code continuously developed during 2010-2014
#
# Changes:
#  MH@2014-02-05
#   changed tcltk progress bar to text based one and modified code,
#   use S4 sts object (no wrapping wanted) and changed to new INLA
#   function name for calculating the transformed marginal.
######################################################################


boda <- function(sts, control=list(range=NULL, X=NULL, trend=FALSE, season=FALSE, prior=c('iid','rw1','rw2'), alpha=0.05, mc.munu=100, mc.y=10, verbose=FALSE,multicore=TRUE)) {

  #Check if the INLA package is available.
  if (!require("INLA")) {
      stop("The boda function requires the INLA package to be installed. The package is not available on CRAN, but can be downloaded by calling source(\"http://www.math.ntnu.no/inla/givemeINLA.R\") as described at http://www.r-inla.org/download in detail.")
  }

  #Possibly speed up the computations by using multiple cores.
  if (is.null(control[["multicore",exact=TRUE]])) { control$multicore <- TRUE }
  if (control$multicore) {
      noCores <- if (require("parallel")) detectCores(logical = TRUE) else 1
      inla.setOption("num.threads", noCores)
  }
  
  #Stop if the sts object is multivariate
  if (ncol(sts)>1) {
      stop("boda currently only handles univariate sts objects.")
  }
  
  # extract data
  observed <- as.vector(observed(sts))
  state <- as.vector(sts@state)
  time <- 1:length(observed)

  
  # clean model data from given outbreaks -- this is now part of the modelling
  #  observed[which(state==1)] <- NA

  ### define range
  # missing range
  if(is.null(control[["range",exact=TRUE]])){ 
    warning('No range given. Range is defined as time from second period until end of time series.')
    control$range <- (sts@freq+1):length(observed)
  }
  # check that range is subset of time series indices
  if(!all(control$range %in% time)){
    stop("Evaluation period 'range' has to be vector of time series indices.")
  }
  #set order of range
  control$range <- sort(control$range) 

  ### show extra output from INLA
  if(is.null(control[["verbose",exact=TRUE]])) {
      control$verbose <- FALSE
  }
  
  ### setting for different models
  if(is.null(control[["trend",exact=TRUE]])){ control$trend <- FALSE }
  if(is.null(control[["season",exact=TRUE]])){ control$season <- FALSE }
  if(!is.logical(control$trend)||!is.logical(control$season)){ 
    stop('trend and season are logical parameters.')
  }
  ### Prior
  prior <- match.arg(control$prior, c('iid','rw1','rw2'))
  if(is.vector(control$X)){
    control$X <- as.matrix(control$X,ncol=1)
  }

  # setting for threshold calcuation
  if(is.null(control[["alpha",exact=TRUE]])){ control$alpha <- 0.05 }
  if(control$alpha <= 0 | control$alpha >= 1){
    stop("The significance level 'alpha' has to be a probability, and thus has to be between 0 and 1.")
  }
  # setting for monte carlo integration
  if(is.null(control[["mc.munu",exact=TRUE]])){ control$mc.munu <- 100 }
  if(is.null(control[["mc.y",exact=TRUE]])){ control$mc.y <- 10 }
  if(!control$mc.munu>0 || control$mc.munu!=round(control$mc.munu,0) || !control$mc.y>0 || control$mc.y!=round(control$mc.y,0)){
    stop('Number of Monte Carlo trials has to be an integer larger than zero')
  }
  
  ### set model formula and data
  modelformula <- paste("observed ~ f(time, model='",prior,"', cyclic=FALSE)", sep="")
  dat <- data.frame(observed=observed, time=time)
  # outbreak id
  if(sum(state)>0){
    modelformula <- paste(modelformula, "+ f(state, model='linear')", sep="")
    dat <- data.frame(dat, state=state) 
  }
  # trend
  if(control$trend){
    modelformula <- paste(modelformula, "+ f(timeT, model='linear')", sep="")
    dat <- data.frame(dat, timeT=time) 
  }
  # season
  if(control$season){
    modelformula <- paste(modelformula, "+ f(timeS, model='seasonal', season.length=",sts@freq,")", sep="")
    dat <- data.frame(dat, timeS=time)
  }
  # covariables
  X.formula <- NULL
  if(!is.null(control$X)){
    if(nrow(control$X)!=length(observed)){
      stop("Argument for covariates 'X' has to have the same length like the time series")
    }
    for(i in 1:ncol(control$X)){
      X.formula <- (paste(X.formula ,'+', colnames(control$X)[i]))
    }
    modelformula <- paste(modelformula, X.formula, sep="")
    dat <- data.frame(dat, control$X)
  }
  modelformula <- as.formula(modelformula)

##### sequential steps #####

  #If there is more than one time point in range, then setup a progress bar 
  #(now text based. Alternative: tcltk based)
  useProgressBar <- length(control$range)>1
  if (useProgressBar) {
    pb <- txtProgressBar(min=min(control$range), max=max(control$range), initial=0,style=3)
  }

  #Allocate vector of thresholds
  xi <- rep(NA,length(observed))

  #Loop over all time points in 'range'
  for(i in control$range){
    # prepare data frame
    dati <- dat[1:i,] 
    dati$observed[i] <- NA #current value to be predicted
    dati$state[i] <- 0 #current state to be predicted

    # fit model and calculate quantile using INLA & MC sampling
#    browser()
    xi[i] <- bodaFit(dat=dati, modelformula=modelformula, prior=prior, alpha=control$alpha, mc.munu=control$mc.munu, mc.y=control$mc.y)

    # update progress bar
    if (useProgressBar) setTxtProgressBar(pb, i)
  }
  # close progress bar
  if (useProgressBar) close(pb)

  # compare observed with threshold an trigger alarm: FALSE=no alarm
  sts@alarm[,1] <- observed > xi 
  sts@upperbound[,1] <- xi
  control$name <- paste('boda(prior=',prior,')',sep='')
  sts@control <- control
  
  # return result as an sts object
  return(sts[control$range,])
}


#######################################################################
# Helper function for fitting the Bayesian GAM using INLA and computing
# the (1-alpha)*100% quantile for the posterior predictive of y[T1]
#
# Parameters:
#  dat - data.frame containing the data
#  modelformula - formula to use for fitting the model with inla
#  prior - what type of prior for the spline c('iid','rw1','rw2')
#  alpha - quantile to compute in the predictive posterior
#  mc.munu - no. of Monte Carlo samples for the mu/size param in the NegBin
#  mc.y - no. of samples for y.
#
# Returns:
#  (1-alpha)*100% quantile for the posterior predictive of y[T1]
######################################################################

bodaFit <- function(dat=dat, modelformula=modelformula,prior=prior,alpha=alpha, mc.munu=mc.munu, mc.y=mc.y,...) {
  
  # set time point
  T1 <- nrow(dat)

  ### fit model
  E <- mean(dat$observed, na.rm=TRUE)
  model <- inla(modelformula,data=dat,family='nbinomial',E=E,control.predictor=list(compute=TRUE))#,verbose=TRUE)
  if(is.null(model)){
    return(qi=NA)
  }
  
  
  # draw sample from marginal posteriori of muT1 & etaT1 to determine predictive
  # quantile by sampling. hoehle: inla.marginal.transform does not exist anymore!  
  # Since the observation corresponding to T1 is NA we manually need to transform
  # the fitted values (had there been an observation this is not necessary!!)
  marg <- try(inla.tmarginal(function(x) exp(x),model$marginals.fitted.values[[T1]]), silent=TRUE)
  #browser()
 
  if(inherits(marg,'try-error')){
      return(qi=NA)
  }
  mT1 <- try(inla.rmarginal(n=mc.munu,marg), silent=TRUE)
  if(inherits(mT1,'try-error')){
      return(qi=NA)
  }
  # take variation in size hyperprior into account by also sampling from it
  mtheta <- model$internal.marginals.hyperpar[[1]]
  theta <- inla.rmarginal(n=mc.munu,mtheta)
  if(inherits(theta,'try-error')){
      return(qi=NA)
  }
  
  #Draw (mc.munu \times mc.y) responses. Would be nice, if we could
  #determine the quantile of the predictive posterior in more direct form
  yT1 <- numeric(mc.munu*mc.y) #NULL
  idx <- seq(mc.y)
  for(j in seq(mc.munu)) {
      idx <- idx + mc.y 
      yT1[idx] <- rnbinom(n=mc.y,size=exp(theta[j]),mu=E*mT1[j])
  }
  
  #Determine the upper (1-alpha)*100% quantile of the predictive posterior
  qi <- quantile(yT1, probs=(1-alpha), type=3, na.rm=TRUE) 
  
  return(qi)
} 

#done bodaFit
