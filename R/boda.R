######################################################################
# An implementation of the Bayesian Outbreak Detection Algorithm (BODA)
# described in Manitz and H{\"o}hle (2013), Biometrical Journal.
#
# Note: The algorithm requires the non-CRAN package INLA to run.
# You can easily install this package as described at
# https://www.r-inla.org/download-install
#
#
# Author:
# The initial code was written by J. Manitz, which was then later
# adapted and modified for integration into the package by M. Hoehle.
# Contributions by M. Salmon and S. Meyer.
#
# Date:
#  Code continuously developed during 2010-2014
#
# Changes:
#  MS@2015-02-18
#   fixed problem that the posterior was drawn from the respective marginals
#   instead of the joint distribution.
#  MH@2014-02-05
#   changed tcltk progress bar to text based one and modified code,
#   use S4 sts object (no wrapping wanted) and changed to new INLA
#   function name for calculating the transformed marginal.
######################################################################


boda <- function(sts, control=list(range=NULL, X=NULL, trend=FALSE, season=FALSE,
                                   prior=c('iid','rw1','rw2'), alpha=0.05, mc.munu=100, 
                                   mc.y=10, verbose=FALSE,
                                   samplingMethod=c('joint','marginals'),
                                   quantileMethod=c("MC","MM"))) {

  #Check if the INLA package is available.
  if (!requireNamespace("INLA", quietly = TRUE)) {
      stop("The boda function requires the INLA package to be installed.\n",
           "  The package is not available on CRAN, but can be easily obtained\n",
           "  from <https://www.r-inla.org/download-install>.")
  }

  #Stop if the sts object is multivariate
  if (ncol(sts)>1) {
      stop("boda currently only handles univariate sts objects.")
  }
  
  # quantileMethod parameter
  if(is.null(control[["quantileMethod",exact=TRUE]])){ 
    control$quantileMethod <- "MC" }
  else {
    control$quantileMethod <- match.arg(control$quantileMethod, 
                                        c("MC","MM"))
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
 # sampling method for the parameters
 samplingMethod <- match.arg(control$samplingMethod, c('joint','marginals'))

  # setting for threshold calculation
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
    pb <- txtProgressBar(min=min(control$range), max=max(control$range), initial=0,
                         style=if (interactive()) 3 else 1)
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
    xi[i] <- bodaFit(dat=dati, samplingMethod=samplingMethod,
                     modelformula=modelformula,
                     alpha=control$alpha, mc.munu=control$mc.munu,
                     mc.y=control$mc.y,
                     quantileMethod=control$quantileMethod,
                     verbose=control$verbose)

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

bodaFit <- function(dat, modelformula, alpha, mc.munu, mc.y,
                    samplingMethod, quantileMethod, verbose = FALSE) {
  
  # set time point
  T1 <- nrow(dat)

  # workaround scoping issue with 'E' in recent versions of INLA
  environment(modelformula) <- environment()

  ### fit model
  link <- 1
  E <- mean(dat$observed, na.rm=TRUE)  # FIXME: is this really needed?
  model <- INLA::inla(modelformula,
                      data=dat,
                      family='nbinomial', E=E, verbose=verbose,
                      control.predictor=list(compute=TRUE,link=link),
                      control.compute=c(list(cpo=FALSE,config=TRUE),
                                        if (packageVersion("INLA") >= "21.07.10")
                                            list(return.marginals.predictor=TRUE)),
                      control.inla = list(int.strategy = "grid",dz=1,diff.logdens = 10))
  if(is.null(model)){ # probably no longer happens in recent versions of INLA
    warning("NULL result from INLA at t = ", T1)
    return(NA_real_)
  }
  
  if(samplingMethod=='marginals'){
    # draw sample from marginal posteriori of muT1 & etaT1 to determine predictive
    # quantile by sampling.
    marg <- model$marginals.fitted.values[[T1]]
    mT1 <- try(INLA::inla.rmarginal(n=mc.munu,marg), silent=TRUE)
    if(inherits(mT1,'try-error')){
        warning("degenerate marginal posterior at t = ", T1)
        return(NA_real_)
    }
    # take variation in size hyperprior into account by also sampling from it
    mtheta <- model$internal.marginals.hyperpar[[1]]
    theta <- exp(INLA::inla.rmarginal(n=mc.munu,mtheta))
  }
  
  if (samplingMethod=='joint'){
    # Sample from the posterior
    ## CAVE: 'model' is not reproducible if num.threads != "1:1" (INLA 22.05.07),
    ##       so there is no point in making the sampling step reproducible
    ##inla.seed <- as.integer(runif(1) * .Machine$integer.max)
    jointSample <- INLA::inla.posterior.sample(
      n = mc.munu, result = model, intern = TRUE,  # seed = inla.seed
      skew.corr = FALSE)  # added with default TRUE in INLA 19.10.30, needs sn
    # take variation in size hyperprior into account by also sampling from it
    theta <- exp(t(sapply(jointSample, function(x) x$hyperpar[[1]])))
    mT1 <- exp(t(sapply(jointSample, function(x) x$latent[[T1]])))
  }
  
  valid <- mT1 >= 0 & theta > 0
  if (any(!valid)) {
    ## a range of (-4.7e-55, 5.8e-52) was seen for mT1 from inla.rmarginal()
    ## which produced an error (-> NA) in previous versions of INLA
    warning("degenerate posterior sampling at t = ", T1)
    return(NA_real_)
    ## mT1 <- mT1[valid]
    ## theta <- theta[valid]
  }
  
  if(quantileMethod=="MC"){
    #Draw (mc.munu \times mc.y) responses. Would be nice, if we could
    #determine the quantile of the predictive posterior in more direct form
    yT1 <- unlist(mapply(rnbinom, size = theta, mu = E*mT1,
                         MoreArgs = list(n = mc.y), SIMPLIFY = FALSE))
    
    qi <- quantile(yT1, probs=(1-alpha), type=3, na.rm=TRUE)
  }
  if(quantileMethod=="MM"){
    minBracket <- qnbinom(p=(1-alpha), 
                          mu=E*min(mT1),
                          size=max(theta))
    
    maxBracket <- qnbinom(p=(1-alpha), 
                          mu=E*max(mT1),
                          size=min(theta))
    
    qi <- qmix(p=(1-alpha), mu=E*mT1, size=theta,
               bracket=c(minBracket, maxBracket))
  }
  return(qi)
} 

#done bodaFit
