###################################################
### chunk number 1: 
###################################################
algo.hhh <- function(disProgObj, control=list(linear=FALSE, nseason=0, period=52, neighbours=FALSE, negbin=FALSE, lambda=TRUE), thetastart=NULL, verbose=TRUE){

  #set default values (if not provided in control)
  if(is.null(control$linear))
    control$linear <- FALSE
  if(is.null(control$nseason))
    control$nseason <- 0
  if(is.null(control$period))
    control$period <- disProgObj$freq #52
  if(is.null(control$neighbours))
    control$neighbours <- FALSE
  if(is.null(control$negbin))
    control$negbin <- FALSE
  if(is.null(control$lambda))
    control$lambda <- TRUE
  
  n <- nrow(disProgObj$observed)
  nareas <- ncol(disProgObj$observed)
  
  #univariate 
  if(nareas ==1){
    control$neighbours <- FALSE
  }
  
  #make "design" matrices
  designRes<- make.design(disProgObj=disProgObj, control=control)
  
  dimtheta <- designRes$dimtheta
  dimbeta <- 2*control$nseason + control$linear   #(gamma_i+delta_i)+beta
  psiIndex <- designRes$psiIndex 
  
  #starting values for optim (without alpha_i's)
  if(!is.null(thetastart)){
    #check dimension of thetastart
    if(length(thetastart) != (dimtheta-nareas)){
      cat('thetastart must be of dimension',dimtheta-nareas,'\n')
     return(NULL)
    }
    theta  <- thetastart
  } else {
    #set starting values for theta   
    #lambda = log(0.5), phi = log(0.1), beta = gamma = delta = 0, psi = 1
    theta <- c(rep(log(0.5),control$lambda), rep(log(0.1),control$neighbours), rep(0, dimbeta), rep(1, control$negbin) )
  }
   
  #starting values for intercepts
  areastart <- log(apply(designRes$Y, 2, sum)/designRes$populationFrac[1,]/(n-1))
  theta <- c(theta, areastart)
  
  # maximize loglikelihood
  mycontrol <- list(fnscale=-1, type=3, maxit=1000)
  suppressWarnings(myoptim <- optim(theta, loglikelihood, control=mycontrol, method="BFGS", hessian=TRUE, designRes=designRes))

  if(myoptim$convergence==0){
    convergence <- TRUE
  } else {
      if(verbose)
        cat("Algorithm has NOT converged. \n")
      res <- list(convergence=FALSE)
      class(res) <- "ah"
      return(res)
  }
  
  loglikelihood <- myoptim$value
  thetahat <- myoptim$par
  fisher <- -myoptim$hessian

  #Added by M.Paul.
  fitted <- meanResponse(thetahat,designRes)     

  #cat("loglik\n")
  #print(loglikelihood)
  #cat("theta\n")
  #print(thetahat)
  #cat("se\n")
  #print(sqrt(-1/diag(myoptim$hessian)))
  #cat("fisher\n")
  #print(fisher,3)
   
  #psi, lambda and phi are on log-scale     
  #-> transformation of estimates, standard errors and fisher (using delta rule)
  #labels for results
  
  # g(theta) = (exp(lambda), exp(phi), beta, gamma, delta, exp(psi), alpha)
  # D is the Jacobian of g
  # D = diag(exp(lambda), exp(phi), 1, 1, 1, exp(psi), 1)
  #
  # mle = g(thetahat)
  # cov = D * fisher^(-1) * D       
  
  thetaNames <- NULL
  D <- NULL
  if(control$lambda){
    thetaNames <- c(thetaNames, "lambda")
    D <- c(D, exp(thetahat[1]))
    thetahat[1] <- exp(thetahat[1])
    if(control$neighbours){
      thetaNames <- c(thetaNames, "phi")
      D <- c(D, exp(thetahat[2]))
      thetahat[2] <-exp(thetahat[2])
    }
  } else if(control$neighbours){ 
    thetaNames <- c(thetaNames, "phi")
    D <- c(D,exp(thetahat[1]))
    thetahat[1] <- exp(thetahat[1])
  }
  
  if(control$linear){
    thetaNames <- c(thetaNames, "beta")
    D <- c(D, 1)
  }
  if(control$nseason > 0){
    gammaDelta <- paste(c("gamma","delta"), rep(1:control$nseason, each=2), sep="")
    thetaNames <- c(thetaNames, gammaDelta )
    D <- c(D, rep(1, length(gammaDelta)))
  }
  if(control$negbin){
    thetaNames <- c(thetaNames, "psi")
    D <- c(D, exp(thetahat[psiIndex]))
    thetahat[psiIndex] <- exp(thetahat[psiIndex])
  }
  
  alpha <- paste("alpha", 1:nareas, sep="")
  thetaNames <- c(thetaNames, alpha)
  D <- c(D, rep(1, length(alpha)))
  
  D <- diag(D)
  
  #fisher <- solve(D)%*%fisher%*%solve(D)   
  
  #Approximation to the inverted fisher info matrix
  cov <- try(D %*% solve(fisher) %*% D, silent=TRUE)
  
  #fisher info is singular
  if(class(cov) == "try-error"){ 
    if(verbose)
      cat('Results are not reliable! Try different starting values. \n')
    res <- list(convergence=FALSE)
    class(res) <- "ah"
    return(res)
  }
  
  if(any(!is.finite(diag(cov))) | any(diag(cov)<0)){   
    if(verbose)        
      cat('Results are not reliable! Try different starting values. \n')
    res <- list(convergence=FALSE)
    class(res) <- "ah"
    return(res)
  }
  
  se <- sqrt(diag(cov))
  
  names(thetahat) <- thetaNames
  names(se) <- thetaNames
  dimnames(cov) <- list(thetaNames,thetaNames)
  
  if(convergence & verbose)
    cat("Algorithm claims to have converged \n")
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #result <- list(thetahat=thetahat, se=se, cov=cov, loglikelihood=loglikelihood, convergence=convergence)
  result <- list(coefficients=thetahat, se=se, cov=cov,
                 loglikelihood=loglikelihood, convergence=convergence,
                 fitted.values=fitted,control=control,disProgObj=disProgObj)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  class(result) <- "ah"
  return(result)
}

print.ah <- function(x,digits = max(3, getOption("digits") - 3), ...){
  if(!x$convergence)
    cat('Results are not reliable! Try different starting values. \n')
  else {
    cat('\nEstimated parameters and standard errors: \n\n')
#~~~~~~~~~~~~~~~~~~~~~~~~
#    print(rbind("Estimates"=x$thetahat, "Std.Error"=x$se),digits=digits)
    print(rbind("Estimates"=coefficients(x), "Std.Error"=x$se),digits=digits)
#~~~~~~~~~~~~~~~~~~~~~~~
    cat('\nloglikelihood:',round(x$loglik,digits=digits-2),'\n\n')
  }
      
}


###################################################
### chunk number 2: 
###################################################
#~~~~~~~~~~~~~~~~~~~~
#################################
# obtain predictions from the fitted algo.hhh model
#
# params:
#  object - a fitted object of class "ah" 
#  newdata - optionally, a disProgObject with which to predict; 
#            if omitted, the fitted mean is returned. 
#  type - the type of prediction required. The default is on the scale of the response 
#         variable (endemic and epidemic part) 
#         the alternative "endemic" returns only the endemic part (i.e. n_it * \nu_it)  
################################
predict.ah<-function(object,newdata=NULL, type=c("response","endemic"),...) {
  type <- match.arg(type,c("response","endemic"))
  control <- object$control

  if(is.null(newdata))
    newdata <- object$disProgObj
  if(!inherits(newdata, "disProg"))
    stop("data must be an object of class disProg\n")

  coefs <- coefficients(object)
  if(type=="endemic"){
    control$lambda <- FALSE
    coefs["lambda"] <- NA
    control$neighbours <- FALSE
    coefs["phi"] <- NA
  }

  design <- make.design(newdata,control=control)

  # in meanResponse the params lambda, phi are "exp()'ed"
  # log() them  to obtain the correct predictions
  coefs["lambda"] <- log(coefs["lambda"])
  coefs["phi"] <- log(coefs["phi"])
  predicted <- meanResponse(coefs[!is.na(coefs)],design)

  return(predicted)
}
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



###################################################
### chunk number 3: 
###################################################
algo.hhh.grid <- function(disProgObj, control=list(linear=FALSE, nseason=0, period=52, neighbours=FALSE, negbin=FALSE, lambda=TRUE), thetastartMatrix, maxTime=1800, verbose=FALSE){

  #set default values (if not provided in control)
  if(is.null(control$linear))
    control$linear <- FALSE
  if(is.null(control$nseason))
    control$nseason <- 0
  if(is.null(control$period))
    control$period <- disProgObj$freq #52
  if(is.null(control$neighbours))
    control$neighbours <- FALSE
  if(is.null(control$negbin))
    control$negbin <- FALSE
  if(is.null(control$lambda))
    control$lambda <- TRUE
  
  nareas <- ncol(disProgObj$observed)
  
  #univariate
  if(nareas ==1){
    control$neighbours <- FALSE
  }
  
  dimthetaStart <- control$lambda+control$neighbours+control$linear+2*control$nseason+control$negbin
  
  if(dimthetaStart == 0){ #only intercepts, grid search not necessary                    
    return(algo.hhh(disProgObj=disProgObj,control=control))
  }
  
  #check dimension of thetastartMatrix
  if(!is.matrix(thetastartMatrix)){
    cat('thetastart must be a matrix with', dimthetaStart, 'columns\n')
    return(NULL)
  }
  if(ncol(thetastartMatrix) != (dimthetaStart)){
    cat('thetastart must be a matrix with',dimthetaStart,'columns\n')
    return(NULL)
  }
  
  #try multiple starting values and return the result with highest likelihood
  #stop search once time limit is exceeded
  i<-0
  nOfIter <- nrow(thetastartMatrix)
  gridUsed <- nOfIter
  
  if(verbose) cat('The size of grid is', nOfIter, '\n')
  
  bestLoglik <- list(loglikelihood = -1e99)
  allLoglik <- matrix(NA,nrow=nOfIter,ncol=1)
  
  while((maxTime > 0) & (i < nOfIter)){
    i <- i+1
    #run algo.hhh with the i-th row of thetastartMatrix as initial values
    time <- system.time(res<-try(algo.hhh(disProgObj=disProgObj,control=control,thetastart=thetastartMatrix[i,],verbose=verbose),silent=!verbose))[3]
    #how much time is left now
    maxTime <- maxTime - time
    
    #print progress information
    if(verbose)
      print(c(niter=i,timeLeft=maxTime,loglik=res$loglikelihood))
    
    #don't consider "useless" results for the search of the best loglikelihood
    if(class(res)!= "try-error" && res$convergence){ 
      #save loglik
      allLoglik[i] <- res$loglikelihood 
      #keep it as bestLoglik if loglikelihood improved  
      if(res$loglikelihood > bestLoglik$loglikelihood){
        bestLoglik <- res
      }
    }
  }
  
  if(maxTime < 0){
    if(verbose) cat('Time limit exceeded, grid search stopped after', i, 'iterations. \n')
    allLoglik <- as.matrix(allLoglik[1:i])
    gridUsed <- i
  }
  
  #algo.hhh did not converge or produced useless results for all starting values, 
  #i.e. there is no result
  if(is.null(coef(bestLoglik))) {
    #convergence <- FALSE
    #cat('Algorithms did not converge, please try different starting values! \n')
    bestLoglik <- list(loglikelihood=NULL,convergence=FALSE)
    #return(NULL)
  }
  
  result <- list(best = bestLoglik, allLoglik = allLoglik,gridSize=nOfIter,gridUsed=gridUsed,convergence=bestLoglik$convergence)
  class(result) <- "ahg"
  return(result) 

}

print.ahg <- function(x, digits = max(3, getOption("digits") - 3),...){
  #how many of the starting values were used
  cat('\nsize of grid: ', x$gridSize, '\n')
  if(x$gridSize!=x$gridUsed)
    cat('grid search stopped after', x$gridUsed, 'iterations \n')

  #algo.hhh did not converge or produced useless results for all starting values,
  #i.e. there is no result
  if(!x$convergence)
    cat('\nAlgorithms did not converge, please try different starting values! \n')
  else {
    #result with highest likelihood
    print.ah(x$best,digits=digits)
  }
}


###################################################
### chunk number 4: 
###################################################
create.grid <- function(params=list(lambda=c(0.1,0.9,5), phi=c(0.1,0.9,5), psi=c(0.3,12,10), beta=c(-0.5,0.5,3), gammaDelta=c(-0.5,0.5,3), nseason=1)){
  
  if(is.null(params$nseason))
    params$nseason <- 1
    
  #number of parameters gamma_i and delta_i
  nOfGammaDelta <- 0
  #negBin model
  negbin <- FALSE
  
  #elements of the grid
  grid <- list()
  
  #autoregressiv parameters lambda, phi and param for NegBin-model psi 
  #must be positive and are used on log-scale in algo.hhh
  if(!is.null(params$psi)){
    grid$psi <- log(seq(params$psi[1],params$psi[2],length=params$psi[3]))
    negbin <- TRUE 
  } 
  if(!is.null(params$lambda))
    grid$lambda <- log(seq(params$lambda[1],params$lambda[2],length=params$lambda[3]))
  if(!is.null(params$phi))
    grid$phi <- log(seq(params$phi[1],params$phi[2],length=params$phi[3]))

  
  #linear trend
  if(!is.null(params$beta))
    grid$beta <- seq(params$beta[1],params$beta[2],length=params$beta[3])
  
  #only one sequence of start values for seasonal parameters gamma AND delta 
  #grid is enlarged later
  if(!is.null(params$gammaDelta)){
    grid$gammaDelta <- seq(params$gammaDelta[1],params$gammaDelta[2],length=params$gammaDelta[3]) 
    #effective number of gamma/delta parameters
    nOfGammaDelta <- 2*params$nseason 
  }
  
  #create grid 
  thetaGrid <- as.matrix(expand.grid(grid))
  gridSize <- nrow(thetaGrid)
  
  #enlargement of grid if seasonal params are present
  if(nOfGammaDelta > 0){
    #one vector of start values for gammaDelta already exists, add it 
    #another (nOfGammaDelta-1)-times to thetaGrid
    gammaDeltaGrid <- matrix(thetaGrid[,"gammaDelta"],nrow=gridSize, ncol=(nOfGammaDelta-1),byrow=FALSE)
    colNames <- paste(c("gamma","delta"), rep(1:(nOfGammaDelta/2), each=2), sep="")
    colnames(gammaDeltaGrid) <- colNames[-1]
    thetaGrid <- cbind(thetaGrid,gammaDeltaGrid)
    colnames(thetaGrid)[colnames(thetaGrid)=="gammaDelta"] <- "gamma1"
  }

  if(negbin){
    #psi is the first column in thetaGrid 
    #in order to get all parameters in the correct order for algo.hhh
    #psi has to be the last one
    psi <- thetaGrid[,"psi"]
    thetaGrid <- thetaGrid[,-1]
    thetaGrid <- cbind(thetaGrid,psi)
  }

  cat('Matrix with starting values for parameters\n',colnames(thetaGrid),'\n')
  return(thetaGrid)
}


###################################################
### chunk number 5: 
###################################################
loglikelihood <- function(theta, designRes){
  control <- designRes$control
  Y <- designRes$Y
  #print(theta,3)                                                            
  mean <- meanResponse(theta=theta, designRes=designRes)
  
  #loglikelihood poisson
  if(control$negbin==FALSE){ 
    result <- sum(dpois(Y, lambda=mean, log=TRUE))  
  }
  
  #loglikelihood negbin
  #psiIndex is position of psi in vector theta
  if(control$negbin==TRUE){ 
    #ensure psi ist positive
    psi <- exp(theta[designRes$psiIndex])

    result <- sum(dnbinom(Y, size=psi, mu=mean, log=TRUE))
  }
  
  return(result)
}


###################################################
### chunk number 6: 
###################################################
meanResponse <- function(theta, designRes){

  control <- designRes$control
  dimtheta <- designRes$dimtheta
  psiIndex <- designRes$psiIndex
 
  Y <- designRes$Y
  X.trendSeason <- designRes$X.trendSeason
  Ym1 <- designRes$Ym1
  Ym1.neighbours <- designRes$Ym1.neighbours
  pop <- designRes$populationFrac

  #check dimension of theta
  if(dimtheta != length(theta)){
    cat('theta must be of dimension',dimtheta,'\n')
    return(NULL)
  }
  
  #unpack parameters and ensure lambda and phi are positive
  lambda <- phi <- 0
  coef.trendSeason <- alpha <- 0  

  #unpack alpha
  alpha <- theta[(psiIndex+1):dimtheta]
  theta <- theta[-((psiIndex+1):dimtheta)]

  #unpack psi
  if(control$negbin){
    psi <- exp(theta[psiIndex])   #psi not needed for mean response, do it anyway
    theta <- theta[-psiIndex]
  }
  #unpack lambda
  if(control$lambda){
    lambda <- exp(theta[1])
    theta <- theta[-1]
  }
  #unpack phi
  if(control$neighbours){
    phi <- exp(theta[1])
    theta <- theta[-1]
  }
  #unpack beta, gamma,delta
  trendSeason <- ((control$nseason >0) | control$linear)
  if(trendSeason){
    coef.trendSeason <- theta
  }
  
  ###################################################################
  # calculation of mean

  #autoregressive part
  #auto=0 if lambda and phi are not used in model
  auto <- Ym1*lambda + Ym1.neighbours*phi

  #trend and seasonal components
  if(trendSeason){
    predtime <- matrix(X.trendSeason%*%coef.trendSeason, byrow=FALSE, ncol=ncol(Y), nrow=nrow(Y))
  } else {
    predtime <- 0
  }
  
  #intercepts for areas
  #matrix with columns (alpha_1,...,alpha_nareas)
  predarea <- matrix(alpha, byrow=TRUE, ncol=ncol(Y), nrow=nrow(Y))

  #results
  mean <- auto + pop*exp(predtime+predarea)

  #Done
  return(mean)
}


###################################################
### chunk number 7: 
###################################################
make.design<-function(disProgObj, control=list(linear=FALSE,nseason=0,period=52,neighbours=FALSE,negbin=FALSE,lambda=TRUE)){

  #set default values (if not provided in control)
  if(is.null(control$linear))
    control$linear <- FALSE
  if(is.null(control$nseason))
    control$nseason <- 0
  if(is.null(control$period))
    control$period <- disProgObj$freq #52
  if(is.null(control$neighbours))
    control$neighbours <- FALSE
  if(is.null(control$negbin))
    control$negbin <- FALSE
  if(is.null(control$lambda))
    control$lambda <- TRUE

  data <- disProgObj$observed
  n <- nrow(data)
  nareas <- ncol(data)

  #theta = (lambda, phi, beta,gamma_i,delta_i,..., psi, alpha_i)
  dimtheta <- control$lambda+control$neighbours+control$linear+2*control$nseason+control$negbin+nareas

  # index of psi (necessary, if negbin=TRUE)
  psiIndex <- dimtheta - nareas
  
  ####################################################################
  # arrange response as matrix
  #Y, Ym1, Ym1.neighbours and population are (nOfobs)x(nOfareas) matrices
  #where nOfobs = n-1 and nOfareas is the number of areas/units
  
  Y <- as.matrix(data[2:n,])
  population <- as.matrix(disProgObj$populationFrac[2:n,])     
  Ym1 <- as.matrix(data[1:(n-1),])
  
  X.trendSeason <- numeric(0)
  X.predict <- numeric(0)
  Ym1.neighbours <- matrix(0,nrow=(n-1),ncol=nareas)    

  # now matrix for neighbours
  if(control$neighbours==TRUE){         
    Ym1.neighbours <- sumNeighbours(disProgObj)[-n,]
  }

  ####################################################################
  # now define design matrix (for trend and seasonality) for each time point
  # definition including prediction!
  
  #trend
  t <- c(1:n)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #t <- t - mean(t)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


  # season design matrix
  season <- NULL
  if(control$nseason > 0){
    omega <- 2*pi/control$period
    for(i in c(1:control$nseason))
      season <- cbind(season,sin(omega*t*i),cos(omega*t*i))
  }

  #combine trend and season  
  if(control$nseason >0)
    X.trendSeason <- season
  if(control$linear)
    X.trendSeason <- cbind(t, season)
    
  if((control$nseason >0) | control$linear){
    #one row with (trend,season)
    X.predict <- X.trendSeason[n,]
    X.trendSeason <- as.matrix(X.trendSeason[c(1:(n-1)),])
  }
             
  result <- list("Y"=Y, "Ym1"=Ym1, "Ym1.neighbours"=Ym1.neighbours, "X.trendSeason"=X.trendSeason, "X.predict"=X.predict, "populationFrac"=population, "dimtheta"=dimtheta,"psiIndex"=psiIndex,"control"=control)

  return(result)
}







