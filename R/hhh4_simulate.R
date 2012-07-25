# helper function to simulate 1 realization from a HHH model
simHHH4 <- function(stsObj, # sts object
                    ar,  # lambda-component of matching size
                    ne,  # phi component
                    end, # nu component
                    psi, # 1/overdisp param(s) if specified or numeric(0)
                    nhood, # weights matrix for phi-component from model
                    start # vector with starting counts
                    ){

  # get required info from sts object
  nUnits <- ncol(stsObj)
  nTime <- nrow(stsObj)

  #################################################
  #Help functions
  ################################################
  psi.inv <- 1/psi
  # draws n random numbers from a NB(mean, psi) distribution
  rNB <- function(n, mean, size = psi.inv){
	rnbinom(n, mu = mean, size=size)
  }

  # weighted sum of all neighbours
  # params: x - vector with counts
  #         nhood - adjacency matrix, 0= no neighbour
  # returns a vector with the sum of "neighbouring counts" for all areas
  wsumN <- function (x, nhood) {
	n <- length(x)
	if(any(nhood>0)){
	  res <- sapply(1:n, function(i) sum(x*nhood[,i]))
	} else {
	  res<- rep(0,n)
	}
	return(res)
  }

  ##################################################

  # simulate from Poisson or NegBin model
  # if psi = 0, then choose Poisson model
  if(length(psi)==0 || psi==0){
	rdistr <- rpois
  } else {
	rdistr<-rNB
  }


  # initialize matrices for the mean mu_i,t and the simulated data x_i,t
  # x_i,0 is set to the mean of n_it*\nu_it
  mu <- matrix(0, ncol=nUnits, nrow=nTime)
  x <- matrix(0, ncol=nUnits, nrow=nTime+1)

  # set starting value to mean observed
  x[1,] <- ifelse(is.null(start),ceiling(colMeans(observed(stsObj))), start)

  # if no ar/ne component then simulate from Poisson/NegBin
  if(all(ar + ne == 0)){
	x <- matrix(rdistr(nUnits*(nTime+1), end), ncol=nUnits, byrow = FALSE)
  } else {
	# simulate data
	for(t in 1:nTime){
	  #mu_i,t = lambda*x_i,t-1 +phi*\sum_j~i wji*x_j,t-1 + nu
	  mu[t,] <- ar[t,] *x[t,] + ne[t,]*wsumN(x[t,], nhood) + end[t,]
	  x[t+1,] <- rdistr(nUnits, mu[t,])
	}
  }

  #return simulated data without first time point
  sts <- stsObj
  observed(sts) <- x[-1,, drop=FALSE]

  return(sts)
}

simulate.ah4 <- function(object, # result from a call to hhh4
                         nsim=1, # number of simulated sts-objects
                         seed=NULL, # seed
                         y.start=NULL,  # starting counts for the epidemic components
                         coefs=NULL,  # if not NULL, these values 
                                      # (in same order [and scale] as coef(object, reparamPsi=TRUE))
                                      # are used to simulate from the model specified in object
                         ...){

  control <- object$control  
  data <- object$stsObj
  # get weights matrix for ne component
  nhood <- control$ne$weights  
  # set to zeros if not specified
  if(is.null(nhood)) nhood <- matrix(0,ncol(data),ncol(data))

  if(!is.null(y.start) && length(y.start)!= ncol(data)){
	stop(paste(sQuote("y.start"), "needs to be of length", ncol(data),".\n"))
  }
  
  # get fitted values of the three components
  model <- interpretControl(control, data)
  theta <- coef(object, reparamPsi=TRUE)  #-> computes 1/exp(logpsi)
  # use different parameter values if specified
  if(!is.null(coefs)){
	if(length(coefs) != length(theta)) stop(paste(sQuote("coefs"), "needs to be of length", length(theta),".\n"))
	names(coefs) <- names(theta)
	theta <- coefs
  }

  # unpack 
  term <- model$terms
  
  pars <- splitParams(theta,model)
  fixed <- pars$fixed
  random <- pars$random
  psi <- pars$overdisp
  
  nGroups <- model$nGroups
  comp <- unlist(term["offsetComp",])
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  
  subset <- 1:nrow(data)#model$subset
   
  toMatrix <- function(par, r=model$nTime, c=model$nUnits){
    matrix(par,r,c, byrow=TRUE)
  }
  
  # go through groups of parameters and compute the lin predictor of each component
  computePartMean <- function(component, subset){
  
    pred <- nullMatrix <- toMatrix(0)
    
    if(!any(comp==component)) return(pred[subset,, drop=FALSE])
    
    for(i in (1:nGroups)[comp==component]){
      fe <- fixed[idxFE==i]
      if(term["unitSpecific",i][[1]]){
        fe <- nullMatrix
        which <- term["which",i][[1]]
        fe[,which] <- toMatrix(fixed[idxFE==i],c=sum(which))
      }
      if(term["random",i][[1]]){
        re <- random[idxRE==i]
        "%m%" <- get(term["mult",i][[1]])
        Z.re <- toMatrix(term["Z.intercept",i][[1]] %m% re)
      } else {
        Z.re <- 0
      }
      X <- term["terms",i][[1]]
      pred <- pred +X*fe + Z.re
    }
    mean <- exp(pred)
    
    return(mean[subset,, drop=FALSE])
  } 
  
  ## autoregressive component
  ar.mean <- computePartMean(1, subset=subset)
  
  ## neighbouring component
  ne.mean <- computePartMean(2, subset=subset)
  
  ## endemic component
  end.mean <- computePartMean(3, subset=subset)*model$offset$end
   
  #now simulate
  if(!is.null(seed)) set.seed(seed)
  if(nsim==1) {
	simData <- simHHH4(data, ar.mean, ne.mean, end.mean, psi, nhood, y.start)
  } else {
	simData <- replicate(nsim, simHHH4(data, ar.mean, ne.mean, end.mean, psi, nhood, y.start), simplify=FALSE)
  }

  #Done
  return(simData) 
}
