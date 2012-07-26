########################################################################
## hhh4 is an extended version of algo.hhh for the sts-class
## The function allows the incorporation of random effects and covariates.
########################################################################

# - some function arguments are currently not used (but will eventually)
# - further argument checks necessary to avoid input errors
# - formula formulation is experimental and not yet implemented in full generality 
# - do some profiling...

hhh4 <- function(stsObj, 
   control = list(
               ar = list(f = ~ -1,        # a formula " exp(x'lamba)*y_t-1 " (ToDo: or a matrix " Lambda %*% y_t-1 ")
                         lag = 1,         # autoregression on y_i,t-lag (currently not used)
                         weights = NULL,  # optional weights, only used if model is a contact matrix (currently not used)
                         initial = NULL   # vector with initial values for parameters if pred is a matrix or if pred = ~1 (not really used ATM)
                         ),
               ne = list(f = ~ -1,        # a formula " exp(x'phi) * \sum_j {w_ji * y_j,t-1} "
                         lag = 1,         # autoregression on y_j,t-lag  (currently not used)
                         weights = NULL,  # weights w_ji, if NULL neighbourhood matrix of stsObj is used
                         initial = NULL   # vector with initial values for parameter if pred = ~1 (not really used ATM)
                         ),
               end = list(f = ~ 1,         # a formula " exp(x'nu) * n_it "
                          offset = NULL,   # optional offset n_it
                          initial = NULL   # vector with initial values for parameter if pred = ~1 (not really used ATM)
                          ),
               family = c("Poisson","NegBin1","NegBinM")[1],
               subset = 2:nrow(stsObj),             # typically 2:nrow(obs) if model contains autoregression
               optimizer = list(tech = "nlminb",    # details for used optimizer (nlminb is default, optim may also be used)
                                stop.tol = 1e-5,
                                stop.niter = 100),
               verbose = FALSE,
               start=list(fixed=NULL,random=NULL,sd.corr=NULL), # list with initials, overrides any initial values in formulas
               data=data.frame(t=epoch(stsObj)-1),    # data.frame or named list with covariates that are specified in the formulas for the 3 components
               keep.terms=FALSE
               )
   ){
        
  #Convert sts objects
  if(inherits(stsObj, "disProg")) stsObj <- disProg2sts(stsObj)
  
  if(is.null(control$data)){
    control$data <- list(.sts=stsObj, t=epoch(stsObj)-1)
  } else {  # coerce to list
    control$data <- modifyList(list(.sts=stsObj),as.list(control$data))
  }

  #set default values (if not provided in control)
  control <- setControl(control, stsObj)
  verbose <- control$verbose
  #* if univariate, no neighbouring term
  #* check nhood matrix
  
  # get model terms
  model <- interpretControl(control,stsObj)
  
  dimFixedEffects <- model$nFE + model$nOverdisp
  dimRandomEffects <- model$nRE
  dimSigma <- model$nSigma


  ##starting values 
  #* -> better default values possible
  theta.start <- model$initialTheta
  Sigma.start <- model$initialSigma
  
  # check if initial values are valid
  # there might be NA's in mu if there are missing values in disProgObj$observed
  mu <- meanHHH(theta.start,model)$mean
  if(any(mu==0, na.rm=TRUE) | any(!is.finite(mu) & !is.na(mu)))
    stop("invalid initial values\n")

  ## maximize loglikelihood 
  cntrl.update <- list(scoreTol=1e-5, paramTol=1e-7,F.inc=0.01, stepFrac=0.5,niter=20) 
  cntrl.stop <- list(tol=control$optimizer$stop.tol,niter=control$optimizer$stop.niter)
  suppressWarnings(myoptim <- fitHHH(theta=theta.start,sd.corr=Sigma.start, model=model, 
                    control=cntrl.stop, cntrl.update=cntrl.update, verbose=verbose, method=control$optimizer$tech))
                 
                 
  if(myoptim$convergence==0){
    convergence <- TRUE
  } else {
      if(verbose)
        cat("Algorithm has NOT converged. \n")
      res <- myoptim
      res$convergence <- FALSE
      class(res) <- "ah4"
      return(res)
  }
  
  loglik <- myoptim$loglik
  
  if(loglik==0){
    if(verbose){
      cat('loglikelihood = 0\n')
      cat('Results are not reliable! Try different starting values. \n')
    }
    res <- myoptim
    res$convergence <- FALSE
    class(res) <- "ah4"
    return(res)
  }

  thetahat <- c(myoptim$fixef,myoptim$ranef)
  fisher <- myoptim$fisher

  # fitted values
  fitted <- meanHHH(thetahat,model)$mean     
 
  if(dimRandomEffects>0){
    Sigma <- myoptim$sd.corr
    Sigma.cov <- solve(myoptim$fisherVar)
    dimnames(Sigma.cov) <- list(names(Sigma),names(Sigma))
  } else {
    Sigma <- Sigma.cov <- NULL
  }

  #Approximation to the inverted fisher info matrix
  cov <- try(solve(fisher), silent=TRUE)

  #fisher info is singular
  if(class(cov) == "try-error"){
    if(verbose){
      cat("Fisher info singular \t loglik=",loglik," \n")
      cat("theta",round(thetahat,2),"\n")
      cat('Results are not reliable! Try different starting values. \n')
    }
    res <- myoptim
    res$convergence <- FALSE
    class(res) <- "ah4"
    return(res)
  }

  if(any(!is.finite(diag(cov))) | any(diag(cov)<0)){
    if(verbose){
      cat("infinite or negative cov\t loglik=",loglik,"\n")
      cat("theta",round(thetahat,2),"\n")
      cat('Results are not reliable! Try different starting values. \n')
    }
    res <- myoptim
    res$convergence <- FALSE
    class(res) <- "ah4"
    return(res)
  }

  se <- sqrt(diag(cov))

  if(convergence & verbose)
    cat("Algorithm converged \n")
  
  margll <- c(marLogLik(Sigma, thetahat,model))
  Sigma.trans <- getSigmai(head(Sigma,model$nVar),tail(Sigma,model$nCorr),model$nVar)
  
  if(control$keep.terms){
	term <- model
  } else {
	term <- NULL
  }

  result <- list(coefficients=thetahat, se=se, cov=cov, 
                 Sigma=Sigma.trans,   # estimated covariance matrix
                 Sigma.orig=Sigma,    # variance parameters on original scale
                 Sigma.cov=Sigma.cov,
                 call=match.call(),
                 dim=c(fixed=dimFixedEffects,random=dimRandomEffects),
                 loglikelihood=loglik, margll=margll, 
                 convergence=convergence,
                 fitted.values=fitted, control=control,terms=term, stsObj=stsObj, 
                 lag=1, nObs=sum(!model$isNA[control$subset,]),nTime=length(model$subset),nUnit=ncol(stsObj))
  
  class(result) <- "ah4"
  return(result)
}

## set default values for model specifications in control
setControl <- function(control, stsObj){
  
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)
  
  
  #set default values (if not provided in control)
  if(is.null( control[["ar", exact = TRUE]] )) {
    control$ar <- list(f = ~ -1, lag = 1, weights = NULL, initial = NULL, isMatrix=FALSE, inModel=FALSE)
  } else {
    if(!is.list(control$ar)){
      stop("\'ar\' must be a list\n")
    } else if(is.null( control$ar[["f", exact = TRUE]] )){
      stop("Predictor for the autregressive component \'ar\' not specified.\n")
    } else if( !( inherits(control$ar$f, "formula") | is.matrix(control$ar$f) ) ){
      stop("\'ar$f\' must be either a formula or a matrix\n")
    } else if( is.matrix(control$ar$f) & any(dim(control$ar$f) != nUnits)){
      stop("\'ar$f\' must be a matrix of size ",nUnits,"\n")  
    }
    if(is.null(control$ar[["lag", exact = TRUE]])) control$ar$lag <- 1
    if(is.null(control$ar[["weights", exact = TRUE]])){
      control$ar$weights <- NULL
    } else if(!is.matrix(control$ar$weights)){
      stop("\'ar$weights\' must be a matrix of size ",nUnits,"\n")          
    } else if(any(dim(control$ar$weights) != nUnits)){
      stop("\'ar$weights\' must be a matrix of size ",nUnits,"\n")    
    }
    if(is.null(control$ar[["initial", exact = TRUE]])) control$ar$initial <- NULL 
    
    # currently only lag=1 considered
    if(control$ar$lag !=1) stop("\'ar$lag\' must be 1\n")
    # check wheter f is a matrix or a formula
    if(is.matrix(control$ar$f)){
      #use identity matrix if weight matrix is missing
      if(is.null(control$ar$weights)) control$ar$weights <- diag(1, nUnits)
      
      control$ar$isMatrix <- TRUE
      control$ar$inModel <- TRUE
    } else {
      control$ar$isMatrix <- FALSE
      if(!is.null(control$ar$weights)) {
        warning("argument \'ar$weights\' is not used\n")
        control$ar$weights <- NULL
      }
      #check if formula is valid
      control$ar$inModel <- isInModel(control$ar$f,"\'ar$f\'") 
    }
  }
  
  ##----------------------------------------------------------------------
  if(is.null( control[["ne", exact = TRUE]] )) {
    control$ne <- list(f = ~ -1, lag = 1, weights = NULL, initial = NULL, inModel = FALSE)
  } else {
    if(!is.list(control$ne)){
      stop("\'ne\' must be a list\n")
    } else if(is.null( control$ne[["f", exact = TRUE]] )){
      stop("Predictor for the neighbour component \'ne\' not specified.\n")
    } else if( !inherits(control$ne$f, "formula") ){
      stop("\'ne$f\' must be a formula\n")
    } 
    if(is.null(control$ne[["lag", exact = TRUE]])) control$ne$lag <- 1
    if(is.null(control$ne[["weights", exact = TRUE]])){
      control$ne$weights <- NULL
    } else if(!is.matrix(control$ne$weights)){
      stop("\'ne$weights\' must be a matrix of size ",nUnits,"\n")        
    } else if(any(dim(control$ne$weights) != nUnits)){
      stop("\'ne$weights\' must be a matrix of size ",nUnits,"\n")    
    }
    if(is.null(control$ne[["initial", exact = TRUE]])) control$ne$initial <- NULL 
    
    # currently only lag=1 considered
    if(control$ne$lag !=1) stop("\'ne$lag\' must be 1\n")
    #check if formula is valid
    control$ne$inModel <- isInModel(control$ne$f,"\'ne$f\'") 
    # if ar$f is a matrix and thus includes neighbouring units, control$ne$f must be empty
    if(control$ar$isMatrix & control$ne$inModel){
      stop("\'ne$f\' is not allowed when \'ar$f\' is a matrix\n")
    }
    
    # if no weight is specified, the nhood-matrix of the stsObj is used
    if(is.null(control$ne$weights)){
      w <- neighbourhood(stsObj)
      diag(w) <- 0
      control$ne$weights <- w
    }
  }
  
  ##----------------------------------------------------------------------
  if(is.null( control[["end", exact = TRUE]] )) {
    control$end <- list(f = ~ 1, offset = 1, initial = NULL, inModel=TRUE)
  } else {
    if(!is.list(control$end)){
      stop("\'end\' must be a list\n")
    } else if(is.null( control$end[["f", exact = TRUE]] )){
      stop("Predictor for the endemic component \'end\' not specified.\n")
    } else if( !inherits(control$end$f, "formula") ){
      stop("\'end$f\' must be a formula\n")
    }   
    if(is.null(control$end[["offset", exact = TRUE]])){
      control$end$offset <- 1
    } else if(!identical(control$end$offset,1) && !is.matrix(control$end$offset)){
      stop("\'end$offset\' must be a matrix of size ",nTime,"x",nUnits,"\n")        
    } else if(!identical(control$end$offset,1) && {(ncol(control$end$offset) != nUnits) | (nrow(control$end$offset) != nTime)}){
      stop("\'end$offset\' must be a matrix of size ",nTime,"x",nUnits,"\n")        
    } else if(any(is.na(control$end$offset))){
      stop("\'end$offset\' cannot contain NA values\n")
    }
    if(is.null(control$end[["initial", exact = TRUE]])) control$end$initial <- NULL
     
    #check if formula is valid
    control$end$inModel <- isInModel(control$end$f,"\'end$f\'") 
     
  }
  
  if(is.null(control[["family",exact=TRUE]]))
    control$family <- "Poisson"
        
  if(is.null(control[["subset",exact=TRUE]])){
    if(nTime <= 2) stop("too few observations\n")
    control$subset <- 2:nTime
  }
        
  if(is.null(control[["optimizer", exact=TRUE]])) control$optimizer <- list(tech="nlminb", stop.tol=1e-5,stop.niter=100)
  if(is.null(control$optimizer[["tech", exact=TRUE]])) control$optimizer$tech <- "nlminb"
  if(is.null(control$optimizer[["stop.tol", exact=TRUE]])) control$optimizer$stop.tol <- 1e-5
  if(is.null(control$optimizer[["stop.niter", exact=TRUE]])) control$optimizer$stop.niter <-100

  if(is.null(control[["verbose",exact=TRUE]]))
    control$verbose <- FALSE
    
  if(is.null(control[["start", exact=TRUE]])){
    control$start <- list(fixed=NULL,random=NULL,sd.corr=NULL)
  } else {
    # set default values (if not provided in start)
    if(is.null(control$start[["fixed",exact=TRUE]]))
      control$start$fixed <- NULL
    if(is.null(control$start[["random",exact=TRUE]]))
      control$start$random <- NULL
    if(is.null(control$start[["sd.corr",exact=TRUE]]))
      control$start$sd.corr <- NULL 
  }
  
  if(is.null(control[["keep.terms",exact=TRUE]]))
    control$keep.terms <- FALSE
        
  control$nTime <- nTime
  control$nUnits <- nUnits

  return(control)
}

# check whether or not one of the three components is included in the model
isInModel <- function(formula, name=deparse(substitute(formula))){
  term <- terms.formula(formula)
  if(attr(term,"response") == 1) stop(name," cannot contain a response variable\n")
  if(attr(term,"intercept") == 1) {
    inModel <- TRUE
  } else {
    # first element is always a list
    inModel <- ifelse(length(attr(term, "variables")) > 1, TRUE, FALSE)
  }
  return(inModel)
}

# used to incorporate covariates and unit-specific effects
fe <- function(x,          # covariate 
               which=NULL, # Null= overall, vector with booleans = unit-specific
               initial=NULL # vector of inital values for parameters
               ){
  sts <- get("env",parent.frame(1))$.sts
  nTime <- nrow(sts)
  nUnits <- ncol(sts)
               
  if(!is.numeric(x)){
    stop("Covariate \'",deparse(substitute(x)),"\' is not numeric\n")
  }
  lengthX <- length(x)
  if(lengthX == 1){
    terms <- matrix(x, nTime, nUnits, byrow=FALSE)
    mult <- "*"
  } else if(lengthX == nTime){
    terms <- matrix(x, nTime, nUnits, byrow=FALSE)
    mult <- "*"  
  } else if(lengthX == nTime*nUnits){
    if(!is.matrix(x)){
     stop("Covariate \'",deparse(substitute(x)),"\' is not a matrix\n")    
    }
    # check dimensions of covariate
    if((ncol(x) != nUnits) | (nrow(x) != nTime)){
     stop("Dimension of covariate \'",deparse(substitute(x)),"\' is not suitably specified\n")    
    }
    terms <- x
    mult <- "*"  
  } else {
    stop("Covariate \'",deparse(substitute(x)),"\' is not suitably specified\n")
  }
  
  intercept <- ifelse(all(terms==1), TRUE, FALSE)
  
  # overall or unit-specific effect?
  unitSpecific <- !is.null(which)
  
  # check argument which
  if(unitSpecific && (!is.vector(which) | (length(which) != nUnits) | !is.logical(which))){
    stop("argument which = \'",deparse(substitute(which)), "\' is not correct\n")
  }
  
  if(unitSpecific){
    terms[,!which] <- 0
  }
  
  # get dimension of parameter
  dim.fe <- ifelse(unitSpecific, sum(which), 1)
  
  #* check length of initial values + set default values?
  if(!is.null(initial) & length(initial) !=dim.fe){
    stop("Initial values for \'",deparse(substitute(x)),"\' must be of length ",dim.fe,"\n")
  } else if(is.null(initial)){
    initial <- rep(0,dim.fe)
  }
  
  summ <- ifelse(unitSpecific,"colSums","sum")
    
  name <- deparse(substitute(x))
  if(unitSpecific) name <- paste(name, colnames(sts)[which], sep=".")
    
  result <- list(terms=terms,
                name=name,
                Z.intercept=NULL,
                which=which,
                dim.fe=dim.fe,
                initial.fe=initial,
                dim.re=0,
                dim.var=0,
                initial.var=NULL,
                initial.re=NULL,
                intercept=intercept,
                unitSpecific=unitSpecific,
                random=FALSE,
                corr=FALSE,
                summ=summ,
                mult=mult
                )
  return(result)
}

# random intercepts
ri <- function(type=c("iid","car")[1], 
               corr=c("none","all")[1],
               initial.var=NULL,  # initial value for variance
               initial.re=NULL 
               ){
  x <- 1
  
  sts <- get("env",parent.frame(1))$.sts
  
  type <- match.arg(type, c("iid","car"))
  corr <- match.arg(corr, c("none","all"))
  corr <- switch(corr, 
                  "none"=FALSE,
                  "all"=TRUE)
                  
  
  nTime <- nrow(sts)
  nUnits <- ncol(sts)
  
  if(type=="iid"){
    Z <- 1
    dim.re <- nUnits 
    mult <- "*"
    
  } else if(type=="car"){
    # construct penalty matrix K
    K <- neighbourhood(sts)
    # check neighbourhood matrix (should contain "1" if i~j and "0" otherwise
    if(any(is.na(K))) stop("neighbourhood matrix contains NA\'s")
    if(!all(K %in% c(0,1))) stop("neighbourhood matrix must contain elements 1 for neighbours and 0 otherwise")
    # number of neighbours
    ne <- colSums(K)
    K <- -1*K
    diag(K) <- ne
    
    dimK <- nrow(K)
    
    # check rank of the nhood, only connected neighbourhoods are allowed
    if(qr(K)$rank != dimK-1) stop("neighbourhood matrix contains islands\n")
    # singular-value decomposition of K
    svdK <- svd(K)
    # just use the positive eigenvalues  of K in descending order 
    # for a the factorisation of the penalty matrix K = LL'
    L <- svdK$u[,-dimK] %*% diag(sqrt(svdK$d[-dimK]))            #* only use non-zero eigenvalues
    
    # Z = L(L'L)^-1
    Z <- L%*%solve(t(L)%*%L)
    
    dim.re <- dimK-1
    mult <- "%*%"
  } 
  
  #* check length of initial values + set default values?
  if(!is.null(initial.re) & length(initial.re) !=dim.re){
    stop("Initial values for \'",type,"\' must be of length ",dim.re,"\n")
  } else if(is.null(initial.re)){
    initial.re <- rnorm(dim.re,0,sd=sqrt(0.001))
  }
  if(!is.null(initial.var) & length(initial.var) !=1){
    stop("Initial values for \'",type,"\' must be of length ",1,"\n")
  } else if(is.null(initial.var)){
    initial.var <- -.5
  }
  
  
  
  result <- list(terms=x,
                name=paste("ri(",type,")",sep=""),
                Z.intercept=Z,
                which=NULL,
                dim.fe=1,
                initial.fe=0,
                dim.re=dim.re,
                dim.var=1,
                initial.var=initial.var,
                initial.re=initial.re,
                intercept=TRUE,
                unitSpecific=FALSE,
                random=TRUE,
                corr=corr,
                summ="colSums",
                mult=mult
                )
  return(result)
}

# check specification of formula
checkFormula <- function(f, env, component){
  term <- terms.formula(f, specials=c("fe","ri"))
  # check if there is an overall intercept
  intercept.all <- attr(term, "intercept") == 1
  
  # list with variables in the component
  vars <- as.list(attr(term,"variables"))
  nVars <- length(vars)-1 # minus 1 because first element is "list"
  
  if(!intercept.all & nVars==0) stop("formula",deparse(substitute(f)),"contains no variables\n")
  
  if(intercept.all){
    res <- c(eval(fe(1),envir=env),list(offsetComp=component))
  } else {
    res <- NULL
  }
  
  # find out fixed effects without "fe()" specification
  fe.raw <- grep("fe(*)|ri(*)", attr(term,"term.labels"), invert=TRUE)
  # evaluate covariates
  if(length(fe.raw)>0){
    for(i in fe.raw)
      res <- cbind(res, c(eval(substitute(fe(x), list(x=vars[[i+1]])),envir=env),list(offsetComp=component)))
  }
  
  # fixed effects
  FE <- attr(term, "specials")$fe
  if(!is.null(FE)){
    for(i in FE){
      res <- cbind(res,c(eval(vars[[i+1]], envir=env),list(offsetComp=component)))
    }  
  }
  
  res <- cbind(res) # ensure res has matrix dimensions
  intercept.fe <- sum(unlist(res["intercept",]))==1
  if(sum(unlist(res["intercept",])) > 1) stop("There can only be one intercept in the formula ",deparse(substitute(f)),"\n")
  
  # random effects
  RI <- attr(term, "specials")$ri
  
  intercept.ri <- !is.null(RI)
  if(length(RI)>1){
    stop("There can only be one random intercept in the formula ",deparse(substitute(f)),"\n")
  }
  if(!is.null(RI)){
    for(i in RI){
      res <- cbind(res,c(eval(vars[[i+1]], envir=env),list(offsetComp=component)))
    }      
  }
  
  
  return(res)
}

# interpret and check the specifications of each component
# control must contain all arguments, i.e. setControl was used
interpretControl <- function(control, stsObj){

  # get environment for evaluation of covariates
  env <- control$data
   
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)

  ##########################################################################
  ##  get the model specifications for each of the three components
  ##########################################################################
  ar <- control$ar
  ne <- control$ne
  end <- control$end
  
  # create response
  Y <- observed(stsObj)
  
  Ym1 <- rbind(matrix(NA,control$ar$lag,nUnits),head(Y, nTime-control$ar$lag)) 
  
  if(control$ne$inModel){
    Ym1.ne <- weightedSumNE(stsObj, control$ne$weights)$neighbours
    Ym1.ne <- rbind(matrix(NA,control$ne$lag,nUnits),head(Ym1.ne, nTime-control$ne$lag))
  } else {
    Ym1.ne <- 1
  }
  
  # list with offsets
  offsets <- list(ar=Ym1, ne=Ym1.ne, end=end$offset)
 
  # determine all NA's
  # offset[[i]] is either 1 or a nTime x nUnits matrix
  isNA <- (is.na(Ym1) | is.na(Ym1.ne) | is.na(Y))
 
 
  ## get terms for all components 
  if(ar$isMatrix){        # ar$f is a matrix:
    stop("not yet implemented\n")
  } else if(ar$inModel){  # ar$f is formula
    all.term <- checkFormula(ar$f,env=env, component=1)    
  } else {
    all.term <-NULL
  }
  
  if(ne$inModel){
    all.term <- cbind(all.term, checkFormula(ne$f,env=env, component=2))
  }
  
  if(end$inModel){
    all.term<- cbind(all.term, checkFormula(end$f,env=env, component=3) )
  }
  
  dim.fe <- sum(unlist(all.term["dim.fe",]))
  dim.re <- sum(unlist(all.term["dim.re",]))
  dim.re.group <- unlist(all.term["dim.re",])
  dim.var <- sum(unlist(all.term["dim.var",]))
  dim.corr <- sum(unlist(all.term["corr",]))
  
  if(dim.corr>0){
    if(dim.var!=dim.corr) stop("Use corr=\'all\' or corr=\'none\' ")
    dim.corr <- switch(dim.corr,0,1,3) 
  }
  
  # the vector with dims of the random effects must be equal if they are correlated
  if(length(unique(dim.re.group[dim.re.group>0]))!=1 & dim.corr>0){ 
    stop("Correlated effects must have same penalty\n")
  }
  
  n <- c("ar","ne","end")[unlist(all.term["offsetComp",])]
  names.fe <- names.var <-NULL
  for(i in 1:length(n)){
    names.fe <- c(names.fe,paste(n[i],all.term["name",i][[1]], sep="."))
    if(all.term["random",i][[1]]) names.var <- c(names.var,paste("sd",n[i],all.term["name",i][[1]], sep="."))
  }
  index.fe <- rep(1:ncol(all.term), times=unlist(all.term["dim.fe",]))
  index.re <- rep(1:ncol(all.term), times=unlist(all.term["dim.re",])) 
  
  fam <- match.arg(control$family,c("Poisson","NegBin1","NegBinM"))
  
  # poisson or negbin model
  if(fam=="Poisson"){
    ddistr <- function(y,mu,size){ dpois(y, lambda=mu,log=TRUE) }
    dim.overdisp <- 0
  } else {
    ddistr <- function(y,mu,size){ dnbinom(y, mu=mu, size=size, log=TRUE) }
    dim.overdisp <- ifelse(fam=="NegBin1",1,nUnits)
  }

  initial.fe <- unlist(all.term["initial.fe",])
  initial.fe.overdisp <- c(initial.fe, rep(2, dim.overdisp))
  initial.re <- unlist(all.term["initial.re",])
  initial.sd<- unlist(all.term["initial.var",])
  initial.sd.corr <- c(initial.sd,rep(0,dim.corr))
  
  # check if a vector of initials is supplied
  if(!is.null(control$start$fixed)){
    if(length(control$start$fixed) != dim.fe+dim.overdisp)
    stop("initial values in start$fixed must be of length ",dim.fe+dim.overdisp,"\n")
    
    initial.fe.overdisp <- control$start$fixed
  }
  if(!is.null(control$start$random)){
  if(length(control$start$random) != dim.re)
    stop("initial values in start$random must be of length ",dim.re,"\n")  
    initial.re <- control$start$random
  }
  if(!is.null(control$start$sd.corr)){
  if(length(control$start$sd.corr) != dim.var+dim.corr)
    stop("initial values in start$sd.corr must be of length ",dim.var+dim.corr,"\n")  
    initial.sd.corr <- control$start$sd.corr
  }

  if(dim.overdisp>1){
    names.overdisp <- paste("overdisp", colnames(stsObj), sep=".")
  } else {
    names.overdisp <- rep("overdisp",dim.overdisp)  # dim.overdisp may be 0
  }
  names(initial.fe.overdisp) <- c(names.fe,names.overdisp)
  initial.theta <- c(initial.fe.overdisp,initial.re)
  names(initial.sd.corr) <- c(names.var,head(paste("corr",1:3,sep="."),dim.corr))
  
  result <- list(response = Y,
                 terms = all.term,
                 nTime = nTime,
                 nUnits = nUnits,
                 nFE = dim.fe,
                 nOverdisp = dim.overdisp,
                 nRE = dim.re,
                 rankRE = dim.re.group,
                 nVar = dim.var,
                 nCorr = dim.corr,
                 nSigma = dim.var+dim.corr,
                 nGroups = ncol(all.term),
                 namesFE = names.fe,
                 indexFE = index.fe,
                 indexRE = index.re,
                 initialTheta = initial.theta,
                 initialSigma = initial.sd.corr,
                 offset = offsets,
                 family=ddistr,
                 subset=control$subset,
                 isNA=isNA
                 )
  return(result)
   
}


getSdCorr <- function(x){
  return(x$Sigma.orig)
}

getCov <- function(x){
  Sigma <- x$Sigma
  corr <- cov2cor(Sigma)
  diag(corr) <- diag(Sigma)
  rownames(corr) <- colnames(corr) <- gsub("sd.","",names(x$Sigma.orig))[grep("sd.",names(x$Sigma.orig))]

  return(corr)
}

 
################################
# Calculates the weighted sum of counts of adjacent areas
# weights are specified in neighbourhood-matrix of the disProgObj
# (experimental atm)
# 
# \nu_i,t = \lambda_y_i,t-1 + \phi*\sum_(j~i) [w_ji*y_j,t-1]
#
# disProgObj$neighbourhood can either be a matrix with weights w_ji (in columns)
# or an array (for time varying weights)
#
# if the neighbourhood-matrix has elements 1 if i~j and 0 otherwise
# weightedSumNeighbours() = sumNeighbours()
###########################################
weightedSumNE <- function(stsObj, weights){

  observed <- observed(stsObj)
  nTime<-nrow(observed)
  nUnits<-ncol(observed)
  neighbours <- matrix(nrow=nTime,ncol=nUnits)
  
  nhood <- weights
  
  #check neighbourhood
  if(any(is.na(nhood)))
    stop("No correct neighbourhood matrix given\n")
  
  ## constant neighbourhood (over time)?
  if(length(dim(nhood))==2){
    # ensure only neighouring areas are summed up
    diag(nhood) <- 0
    nhood <- array(nhood,c(nUnits,nUnits,nTime))
  
  } else if(length(dim(nhood))==3){
    if(any(dim(nhood)[1:2]!= nUnits) | dim(nhood)[3] != nTime) 
      stop("neighbourhood info incorrect\n")
  }
    
  # number of neighbours
  nOfNeighbours <-colSums(nhood[,,1]>0)
    
  for(i in 1:ncol(observed)){
    weights <- t(nhood[,i,])
    weightedObs <- observed*weights
    neighbours[,i] <- rowSums(weightedObs, na.rm=TRUE)
  }
  
  return(list(neighbours=neighbours, nOfNeighbours=nOfNeighbours))
}


splitParams <- function(theta, model){
  fixed <- head(theta,model$nFE)
  random <- tail(theta, model$nRE)
  if(model$nOverdisp >0){
    overdisp <- theta[(model$nFE+1):(model$nFE+model$nOverdisp)]
  } else  overdisp <- numeric(0)
  
  return(list(fixed=fixed,random=random,overdisp=overdisp))
}

# compute predictor
meanHHH <- function(theta, model){

  # unpack 
  Y <- model$response
  term <- model$terms
  offsets <- model$offset
  
  pars <- splitParams(theta,model)
  fixed <- pars$fixed
  random <- pars$random
  
  nGroups <- model$nGroups
  comp <- unlist(term["offsetComp",])
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  
  subset <- model$subset
  
  #* check at higher level
  if(length(theta) != model$nRE+model$nFE+model$nOverdisp){
    stop("theta must be of length ",model$nRE+model$nFE+model$nOverdisp,"\n")
  }
  
  toMatrix <- function(par, r=model$nTime, c=model$nUnits){
    matrix(par,r,c, byrow=TRUE)
  }
  
  # set missing values in observed Y to NA
  # offset[[i]] is either 1 or a nTime x nUnits matrix
  isNA <- model$isNA

  # go through groups of parameters and compute the lin predictor of each component
  computePartMean <- function(component, isNA, subset){
  
    pred <- nullMatrix <- toMatrix(0)
    pred[isNA] <- NA
    
    if(!any(comp==component)) return(pred[subset,])
    
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
    mean <- exp(pred)*offsets[[component]]
    
    mean[isNA] <-NA
    return(mean[subset,])
  } 
  
  ## autoregressive component
  ar.mean <- computePartMean(1, isNA=isNA, subset=subset)
  
  ## neighbouring component
  ne.mean <- computePartMean(2, isNA=isNA, subset=subset)
  
  ## endemic component
  end.mean <- computePartMean(3, isNA=isNA, subset=subset)
   
  #results
  mean <- ar.mean + ne.mean + end.mean

  #Done
  return(list(mean=mean,epidemic=ar.mean+ne.mean,endemic=end.mean,epi.own=ar.mean,epi.neighbours=ne.mean))
}

############################################
penLogLik <- function(theta, sd.corr, model){

  # unpack 
  Y <- model$response[model$subset,,drop=FALSE]
  ddistr <- model$family
  
  sd <- head(sd.corr,model$nVar)
  corr <- tail(sd.corr,model$nCorr)

  
  pars <- splitParams(theta,model)
  randomEffects <- pars$random
  
  #ensure overdispersion param is positive
  overdispParam <- exp(pars$overdisp)
    
  if(model$nOverdisp > 1) {
    overdispParam <- matrix(overdispParam,ncol=model$nUnits, nrow=model$nTime, byrow=TRUE)[model$subset,,drop=FALSE]
  }
  
  mu <- meanHHH(theta=theta, model=model)$mean
  
  #######    
  lpen <- 0
  # check if there are random effects
  if(model$nVar >0){
    dimBlock<- model$rankRE[model$rankRE>0]
    Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
    
    lpen <- -.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
    
  } 
  
  #loglikelihood 
  ll.units <- colSums(ddistr(Y,mu,overdispParam), na.rm=TRUE)
  
  ll <- sum(ll.units) + lpen
  attr(ll, "loglik") <- ll.units
  attr(ll, "logpen") <- lpen
  return(ll)
}


penScore <- function(theta, sd.corr, model){
  
  if(any(is.na(theta) | !is.finite(theta))){ 
    return(rep(NA,length(theta)))
  }
  
  # unpack 
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  nUnits <- model$nUnits
  
  term <- model$terms
  offsets <- model$offset
  
  nGroups <- model$nGroups
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  
  sd <- head(sd.corr,model$nVar)
  corr <- tail(sd.corr,model$nCorr)
    
  dimFE <- model$nFE
  dimRE <- model$nRE
  
  # predictor
  mu <- meanHHH(theta=theta, model=model)
  mu.end <- mu$endemic
  mu.ar <- mu$epi.own
  mu.ne <- mu$epi.neighbours
  meanTotal <- mu$mean

  pars <- splitParams(theta,model)
  #ensure overdispersion param is positive
  psi <- exp(pars$overdisp)
  if(model$nOverdisp > 1) {
    psi <- matrix(psi,ncol=model$nUnits, nrow=model$nTime, byrow=TRUE)[subset,,drop=FALSE]
  }
  #random effects
  randomEffects <- pars$random
  dimBlock<- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)

  ############################################################
  ## helper function for derivatives:
  # negbin model or poisson model
  if(model$nOverdisp > 0){
    psiPlusMu <- psi + meanTotal
    
    # helper function for derivatives: negbin
    derivHHH <- function(dmu){
      (-psi/psiPlusMu +Y/meanTotal -Y/psiPlusMu)*dmu
    }
      
  } else {
    # helper function for derivatives: poisson
    derivHHH <- function(dmu){
      Y *(dmu/meanTotal) - dmu
    }

  }

  ############################################################

  # go through groups of parameters and compute the gradient of each component
  computeGrad <- function(mean.comp, subset){
  
    grad.fe <- NULL
    grad.re <- NULL
       
    for(i in 1:nGroups){
      comp <- term["offsetComp",i][[1]]
      Xit<- term["terms",i][[1]] # eiter 1 or a matrix with values
      if(is.matrix(Xit)){
        Xit <- Xit[subset,,drop=FALSE]
      }
      summ <- get(term["summ",i][[1]])
      dTheta <- derivHHH(mean.comp[[comp]]*Xit)
      dTheta[isNA] <- 0
      
      if(term["unitSpecific",i][[1]]){
        which <- term["which",i][[1]]
        dTheta <- summ(dTheta,na.rm=TRUE)[ which ]
        grad.fe <- c(grad.fe,dTheta)
        
      } else if(term["random",i][[1]]){
        Z <- term["Z.intercept",i][[1]]  
        "%m%" <- get(term["mult",i][[1]])
        dRTheta <- colSums(dTheta %m% Z)  #dTheta must not contain NA's (set NA's to 0)      
        grad.re <- c(grad.re,dRTheta)
        grad.fe <- c(grad.fe, sum(dTheta))
      } else{
        grad.fe <- c(grad.fe,summ(dTheta,na.rm=TRUE))
      }
      
   }
    
    return(list(fe=grad.fe,re=grad.re))
  } 
  
  gradients <- computeGrad(mean.comp=list(ar=mu.ar, ne=mu.ne,end=mu.end),subset=subset)
  
  
  # gradient for overdispersionParameter psi
  if(model$nOverdisp > 0){
    dPsi <- psi*(psigamma(Y+psi)-psigamma(psi) +log(psi)+1 - log(psiPlusMu) -psi/psiPlusMu -Y/psiPlusMu)
    
    # multiple psi_i's or single psi?
    if(model$nOverdisp > 1){
      grPsi <- colSums(dPsi, na.rm=TRUE)
    }else {
      grPsi <- sum(dPsi, na.rm=TRUE)
    }
    
    if(any(is.na(grPsi))){
      warning("derivatives for psi not computable\n")
      return(rep(NA,length(theta)))
    }
    
  } else {
    grPsi <- NULL  
  }
 
  ##------------------------------------------------
  ## Random Effects
  ##------------------------------------------------
  if(dimRE >0){
    # penalty of score function
    s.pen <- Sigma.inv %*% randomEffects
    # unpenalized score function for random effects
    grRandom <- gradients$re
    
    if(length(grRandom) != length(s.pen))
      cat("length of s(b) != Sigma.inv %*%b \n")
      
    grRandom <- c(grRandom - s.pen)
  } else {
    grRandom <- NULL
  }
  res <- c(gradients$fe, grPsi ,grRandom)
  
  return(res)

}

penFisher <- function(theta, sd.corr, model, attributes=FALSE){
   
  # unpack 
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
    
  nUnits <- model$nUnits
  
  term <- model$terms
  offsets <- model$offset
  
  nGroups <- model$nGroups
  idxFE <- model$indexFE
  idxRE <- model$indexRE
    
  dimFE <- model$nFE
  dimRE <- model$nRE
  dimPsi <- model$nOverdisp
  
  sd <- head(sd.corr,model$nVar)
  corr <- tail(sd.corr,model$nCorr)
  
  # predictor
  mu <- meanHHH(theta=theta, model=model)
  mu.end <- mu$endemic
  mu.ar <- mu$epi.own
  mu.ne <- mu$epi.neighbours
  meanTotal <- mu$mean

  pars <- splitParams(theta,model)
  #ensure overdispersion param is positive
  psi <- exp(pars$overdisp)
  if(model$nOverdisp > 1) {
    psi <- matrix(psi,ncol=model$nUnits, nrow=model$nTime, byrow=TRUE)[subset,,drop=FALSE]
  }
  #random effects
  randomEffects <- pars$random
  dimBlock<- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  
 
  ## helper functions for derivatives:
  # negbin model or poisson model
  if(model$nOverdisp>0){
    psiPlusY <- psi + Y
    psiPlusMu <- psi + meanTotal
    psiPlusMu2 <- psiPlusMu^2
    
    # helper function for derivatives: negbin
    deriv2HHH <- function(dTheta_l,dTheta_k,dTheta_lk){
        dTheta_l*dTheta_k*(psi/psiPlusMu2 - Y/(meanTotal^2) + Y/psiPlusMu2) + dTheta_lk*(-psi/psiPlusMu +Y/meanTotal - Y/psiPlusMu) 
    }
    
    dThetadPsi <- function(dTheta){
        psi*(-dTheta/psiPlusMu + (psi + Y)*dTheta/psiPlusMu2)
    }
    
    dPsidPsi <- function(){
        dPsi <- psi*(psigamma(psiPlusY)-psigamma(psi) +log(psi)+1 - log(psiPlusMu) -psi/psiPlusMu -Y/psiPlusMu)
        psi*(psigamma(psiPlusY,1)*psi - psigamma(psi,1)*psi + 1 - psi/psiPlusMu - psi*(meanTotal-Y)/psiPlusMu2) +dPsi*psi
    }
    
  } else {
    # helper function for derivatives: poisson
    deriv2HHH <- function(dTheta_l,dTheta_k,dTheta_lk){
        -Y*dTheta_l*dTheta_k/(meanTotal^2) + (Y/meanTotal-1)*dTheta_lk
    }
  }
     
  # go through groups of parameters and compute the hessian of each component
  computeFisher <- function(mean.comp, subset){
    # initialize hessian
    hessian.FE.FE <- matrix(0,dimFE,dimFE)
    hessian.FE.RE <- matrix(0,dimFE,dimRE)
    hessian.RE.RE <- matrix(0,dimRE,dimRE)
    
    hessian.FE.Psi <- matrix(0,dimFE, dimPsi)
    hessian.Psi.RE <- matrix(0,dimPsi, dimRE+dimPsi)

    ##
    i.fixed <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]  
        "%m%" <- get(term["mult",j][[1]])
        dIJ <- colSums(didj %m% Z.j)      # didj must not contain NA's (all NA's set to 0)
        hessian.FE.RE[idxFE==i,idxRE==j] <<- dIJ
        dIJ <-sum(dIJ)
        
      } else if(unitSpecific.j){
        dIJ <- colSums(didj,na.rm=TRUE)[ which.j ]
      } else {
        dIJ <- sum(didj,na.rm=TRUE)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ  
    }
    ##
    i.unit <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]] 
        "%m%" <- get(term["mult",j][[1]]) 
        dIJ <- colSums(didj %m% Z.j)     # didj must not contain NA's (all NA's set to 0)
        hessian.FE.RE[idxFE==i,idxRE==j] <<- diag(dIJ)[ which.i, ]
        
      } else if(unitSpecific.j){
        which.ij <- (which.i & which.j)
        dIJ <- diag(colSums(didj,na.rm=TRUE))[ which.i, which.j ] 
      } else {
        dIJ <- colSums(didj,na.rm=TRUE)[ which.i ]
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ  
    }
    ##
    i.random <- function(){
      if(random.j){
        dIJ <- sum(didj)
        Z.j <- term["Z.intercept",j][[1]]  
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums( didj %mj% Z.j)
        
        if(length(Z.j)==1 & length(Z.i)==1){
          Z <- Z.i*Z.j
          hessian.RE.RE[which(idxRE==i),idxRE==j] <<- diag(colSums( didj %m% Z))
        } else if(length(Z.j)==1 & length(Z.i)>1){         #*
          Z.j <- diag(1,model$nUnits)
          for(k in 1:ncol(Z.j)){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[idxRE==i,which(idxRE==j)[k]] <<- colSums( didj %m% Z)
          }  
        } else if(length(Z.j)>1 & length(Z.i)==1){         #* 
          Z.i <- diag(1,model$nUnits)
          for(k in 1:ncol(Z.i)){
            Z <- Z.i[,k]*Z.j
            hessian.RE.RE[which(idxRE==i)[k],idxRE==j] <<- colSums( didj %mj% Z)
          }  
        } else {         
          for(k in 1:ncol(Z.j)){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[which(idxRE==i)[k],idxRE==j] <<- colSums( didj %m% Z)
          }  
        }        
              
      } else if(unitSpecific.j){
        dIJ <-  colSums(didj)[ which.j ]
      } else {
        dIJ <- sum(didj,na.rm=TRUE)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ    
    }
  ##----------------------------------------------           

    for(i in 1:nGroups){ #go through rows of hessian
      # parameter group belongs to which components
      comp.i <- term["offsetComp",i][[1]]      
      # get covariate value
      Xit<- term["terms",i][[1]] # eiter 1 or a matrix with values
      if(is.matrix(Xit)){
        Xit <- Xit[subset,,drop=FALSE]
      }
      m.Xit <- mean.comp[[comp.i]]*Xit
      
      random.i <- term["random",i][[1]]
      unitSpecific.i <- term["unitSpecific",i][[1]]
      
      ## d l(theta,x) /dpsi dpsi
      if(dimPsi==1){
         hessian.Psi.RE[1,1] <- sum(dPsidPsi(), na.rm=TRUE)
      } else if(dimPsi >1){
         hessian.Psi.RE[, 1:dimPsi] <- diag(colSums(dPsidPsi(), na.rm=TRUE))
      } 
      
      if(random.i){
        Z.i <- term["Z.intercept",i][[1]]   
        "%m%" <- get(term["mult",i][[1]])
         
        fillHess <-i.random
        if(dimPsi==1){
          dPsi<- dThetadPsi(m.Xit)
          dPsi[isNA] <- 0
          hessian.FE.Psi[idxFE==i,]<- sum(dPsi)
          hessian.Psi.RE[,c(FALSE,idxRE==i)] <- colSums( dPsi%m%Z.i)
        } else if(dimPsi>1){
          dPsi<- dThetadPsi(m.Xit)
          dPsi[isNA] <- 0
          hessian.FE.Psi[idxFE==i,]<- colSums( dPsi)
          hessian.Psi.RE[,c(FALSE,idxRE==i)] <- diag(colSums( dPsi%m%Z.i))
        }
        
      } else if(unitSpecific.i){
        which.i <- term["which",i][[1]]
        fillHess <- i.unit
        if(dimPsi>0){
          dPsi<- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
          if(dimPsi==1){ 
            hessian.FE.Psi[idxFE==i,] <- dPsi[which.i] 
          } else {
            hessian.FE.Psi[idxFE==i,] <- diag(dPsi)[which.i,]
          }
        }
        
      } else {
        fillHess <- i.fixed
        if(dimPsi>0){
          dPsi<- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
          if(dimPsi==1){ dPsi <- sum(dPsi) }
          hessian.FE.Psi[idxFE==i,] <-dPsi
        }
      }
            

      for(j in 1:nGroups){
        comp.j <- term["offsetComp",j][[1]]
        
        Xjt<- term["terms",j][[1]] # eiter 1 or a matrix with values
        if(is.matrix(Xjt)){
          Xjt <- Xjt[subset,,drop=FALSE]
        }
        # if param i and j do not belong to the same component, d(i)d(j)=0
        if(comp.i!=comp.j){
          m.Xit.Xjt <- 0
        } else {
          m.Xit.Xjt <- m.Xit*Xjt        
        }
        
        didj <- deriv2HHH(dTheta_l=m.Xit,dTheta_k=mean.comp[[comp.j]]*Xjt, dTheta_lk=m.Xit.Xjt)
        didj[isNA]<-0

        random.j <- term["random",j][[1]]      
        unitSpecific.j <- term["unitSpecific",j][[1]]
        which.j <- term["which",j][[1]]
        
        fillHess()
      }        
  
    }
    
    
    ######################################################### 
    ## fill lower triangle of hessians and combine them
    ########################################################
    hessian <- rbind(cbind(hessian.FE.FE,hessian.FE.Psi,hessian.FE.RE),
                    cbind(matrix(0,dimPsi,dimFE),hessian.Psi.RE),
                    cbind(matrix(0,dimRE,dimFE+dimPsi),hessian.RE.RE))
          
    diagHessian <- diag(hessian)
    hessian[lower.tri(hessian)] <- 0    #* only use upper.tri?
    fisher <- -(hessian + t(hessian))
    diag(fisher) <- -diagHessian   
  
    return(fisher)
  } 

  fisher<- computeFisher(mean.comp=list(ar=mu.ar, ne=mu.ne,end=mu.end),subset=subset)
 
  
  if(dimRE >0){
    ## Add penalty for Random Effects
    pen <- matrix(0, dimFE+dimPsi+dimRE, dimFE+dimPsi+dimRE)
    pen[-(1:(dimFE+dimPsi)),-(1:(dimFE+dimPsi))] <- Sigma.inv
    Fpen <- fisher + pen
  } else {
    pen <- 0
    Fpen <- fisher
  }

  if(attributes){
    attr(Fpen, "fisher") <- fisher
    attr(Fpen, "pen") <- pen
    return(Fpen)
  } else return(Fpen)

}

#################################################
sqrtOf1pr2 <- function(r){ sqrt(1+r^2)}

getSigmai <- function(sd,    # vector of length dim with stdev's
                      correlation, # vector of length dim with correlation parameters, =0 if un-correlated
                      dim
                      ){
  if(dim==0) return(NULL)
                      
  D <- diag(exp(sd),dim)
  if(length(correlation)>0){
    L <- diag(1,dim,dim)
    L[2,1:2] <- c(correlation[1],1)/sqrtOf1pr2(correlation[1])
    if(dim==3){
      L[3,] <- c(correlation[2:3],1)/sqrtOf1pr2(correlation[2])
      L[3,2:3] <- L[3,2:3]/sqrtOf1pr2(correlation[3])
    }
  } else{
    L <- diag(1,dim)
  }
  Sigma.i <- D%*%L%*%t(L)%*%D
  return(Sigma.i)                     
}

getSigmaiInv <- function(sd,    # vector of length dim with stdev's
                      correlation, # vector of length dim with correlation parameters, =0 if un-correlated
                      dim
                      ){
                      
  if(dim==0) return(NULL)                     
  r <- correlation
                     
  Dinv <- diag(exp(-sd),dim)
  if(length(correlation)>0){
    L <- diag(1,dim,dim)
    L[2,1:2] <- c(-r[1],sqrtOf1pr2(r[1]))
    if(dim==3){
      L[3,1] <- r[1]*r[3]-r[2]*sqrtOf1pr2(r[3])
      L[3,2] <- -L[2,2]*r[3]
      L[3,3] <- sqrtOf1pr2(r[2])*sqrtOf1pr2(r[3])
    }
  } else{
    L <- diag(1,dim)
  }

  Sigma.i.inv <- Dinv%*%t(L)%*%L%*%Dinv   
  return(Sigma.i.inv)                     
}

#* allow blockdiagonal matrix blockdiag(A,B), with A=kronecker product, B=diagonal matrix?
getSigmaInv <- function(sd, correlation, dimSigma, dimBlocks, SigmaInvi=NULL){
  if(is.null(SigmaInvi)){
    SigmaInvi<- getSigmaiInv(sd,correlation,dimSigma)
  }
  if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    SigmaInv <- kronecker(SigmaInvi,diag(1,dimBlocks[1]))
  } else { # kronecker product not possible -> correlation=0 
    SigmaInv <- diag(rep(diag(SigmaInvi),dimBlocks))
  }
  return(SigmaInv)
}

getSigma <- function(sd, correlation, dimSigma, dimBlocks, Sigmai=NULL){
  if(is.null(Sigmai)){
    Sigmai <- getSigmai(sd,correlation,dimSigma)
  }
  if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    Sigma <- kronecker(Sigmai,diag(1,dimBlocks[1]))
  } else { # kronecker product not possible -> correlation=0 
    Sigma <- diag(rep(diag(Sigmai),dimBlocks))
  }
  return(Sigma)
}

# approximate marginal likelihood for variance components
marLogLik <- function(sd.corr, theta,  model, fisher.unpen=NULL){
  
  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(Inf)
  }
  
    sd <- head(sd.corr,dimVar)
    corr <- tail(sd.corr,dimCorr)
  
  if(any(is.na(sd.corr))){
   cat("WARNING: NAs in variance components\n") 
    return(NA)      
  }
  
  
  dimFE <- model$nFE
  dimOver <- model$nOverdisp
  dimFE.O <- dimFE+dimOver
  dimRE <- model$nRE
  
  dimBlocks<- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)
  
  
  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attributes(penFisher(theta, sd.corr, model,attributes=TRUE))$fisher
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen 
  fisher[-(1:dimFE.O),-(1:dimFE.O)] <- fisher[-(1:dimFE.O),-(1:dimFE.O)] + Sigma.inv
  
  F.inv <- try(solve(fisher),silent=TRUE)
  
  if(inherits(F.inv,"try-error")){
    cat("\n WARNING: penalized Fisher is singular!\n")
    return(NA)  
  }
   
  F.inv.RE <- F.inv[-(1:dimFE.O),-(1:dimFE.O)]
  
  pars <- splitParams(theta,model)  
  randomEffects <- pars$random

  
  # penalized part of likelihood
  # compute -0.5*log(|Sigma|) - 0.5*RE' %*% Sigma.inv %*% RE
  # where -0.5*log(|Sigma|) = -dim(RE_i)*[Sum(sd_i) -0.5*log(1+corr_i^2)]
  loglik.pen <- sum(-dimBlocks*sd) - 0.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
  if(dimCorr >0){
    loglik.pen <- loglik.pen + 0.5*dimBlocks[1]*sum(log(1+corr^2))
  }
  
  ## approximate marginal likelihood
  lmarg <- loglik.pen -0.5*c(determinant(fisher,logarithm=TRUE)$modulus)
      
  return(lmarg)
}


marScore <- function(sd.corr, theta,  model, fisher.unpen=NULL){
  
  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(rep(NA,dimVar))
  }
  
    sd <- head(sd.corr,dimVar)
    corr <- tail(sd.corr,dimCorr)
  
  if(any(is.na(sd.corr))){
   cat("WARNING: NAs in variance components\n") 
    return(rep(NA,dimVar))
  }
  
  
  dimFE <- model$nFE
  dimOver <- model$nOverdisp
  dimFE.O <- dimFE+dimOver
  dimRE <- model$nRE
  
  dimBlocks<- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)
  
  
  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attributes(penFisher(theta, sd.corr, model,attributes=TRUE))$fisher
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen 
  fisher[-(1:dimFE.O),-(1:dimFE.O)] <- fisher[-(1:dimFE.O),-(1:dimFE.O)] + Sigma.inv
  
  F.inv <- try(solve(fisher),silent=TRUE)
  
  if(inherits(F.inv,"try-error")){
    cat("\n WARNING: penalized Fisher is singular!\n")
    return(rep(NA,dimVar))
  }
   
  F.inv.RE <- F.inv[-(1:dimFE.O),-(1:dimFE.O)]
  
  pars <- splitParams(theta,model)  
  randomEffects <- pars$random
   
  ## compute marginal score and fisher for each variance component
  # initialize score and fisher info
  marg.score <- rep(NA,dimSigma)
  
  ## specify functions for derivatives
  deriv1 <- switch(dimVar,dSigma1, dSigma2, dSigma3)
  
  d1Sigma <- deriv1(sd, corr)
  Sigmai.inv <- getSigmaiInv(sd, corr, dimVar)
  
 
  # derivation of log determinant
  # -.5*tr(Sigma^-1 %*% dSigma/ds) = -R (for sd.i) 
  #                                =  R*corr.i/(corr.i^2+1) (for corr.i)
  d1logDet <- c(-dimBlocks,dimBlocks[1]*corr/(corr^2+1))
  
  # go through all variance components
  for(i in 1:dimSigma){
    dSi <- -Sigmai.inv %*% d1Sigma[,,i] %*% Sigmai.inv
    dS.i <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSi)
    marg.score[i] <- d1logDet[i] - 
                      0.5* t(randomEffects)%*% dS.i %*% randomEffects -
                      0.5* sum(diag(F.inv.RE %*% dS.i))                     
  }  
  
  return(marg.score)
}

marFisher <- function(sd.corr, theta,  model, fisher.unpen=NULL){
  
  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(matrix(NA,dimVar,dimVar))   
  }
  
    sd <- head(sd.corr,dimVar)
    corr <- tail(sd.corr,dimCorr)
  
  if(any(is.na(sd.corr))){
   cat("WARNING: NAs in variance components\n") 
    return(matrix(NA,dimVar,dimVar))   
  }
  
  
  dimFE <- model$nFE
  dimOver <- model$nOverdisp
  dimFE.O <- dimFE+dimOver
  dimRE <- model$nRE
  
  dimBlocks<- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)
  
  
  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attributes(penFisher(theta, sd.corr, model,attributes=TRUE))$fisher
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen 
  fisher[-(1:dimFE.O),-(1:dimFE.O)] <- fisher[-(1:dimFE.O),-(1:dimFE.O)] + Sigma.inv
  
  F.inv <- try(solve(fisher),silent=TRUE)
  
  if(inherits(F.inv,"try-error")){
    cat("\n WARNING: penalized Fisher is singular!\n")
    return(matrix(NA,dimVar,dimVar))   
  }
   
  F.inv.RE <- F.inv[-(1:dimFE.O),-(1:dimFE.O)]
  
  pars <- splitParams(theta,model)  
  randomEffects <- pars$random

  marg.hesse <- matrix(NA,dimSigma,dimSigma)
  
  ## specify functions for derivatives
  deriv1 <- switch(dimVar,dSigma1, dSigma2, dSigma3)
  deriv2 <- switch(dimVar,d2Sigma1, d2Sigma2, d2Sigma3)
  
  d1Sigma <- deriv1(sd, corr)
  d2Sigma <- deriv2(sd, corr, d1Sigma)
  Sigmai.inv <- getSigmaiInv(sd, corr, dimVar)
  
 
  # 2nd derivatives of log determinant
  d2logDet <- diag(c(rep(0,dimVar),-dimBlocks[1]*(corr^2-1)/(corr^2+1)^2),dimSigma)
  
  # go through all variance components
  for(i in 1:dimSigma){
	# compute first derivative of the penalized Fisher info (-> of Sigma^-1) 
    # with respect to the i-th element of Sigma (= kronecker prod. of Sigmai and identity matrix)
	# Harville Ch15, Eq. 8.15: (d/d i)S^-1 = - S^-1 * (d/d i) S * S^-1
    SigmaiInv.d1i <- Sigmai.inv %*% d1Sigma[,,i]
    dSi <- -SigmaiInv.d1i %*% Sigmai.inv
    dS.i <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSi)
    
	# compute second derivatives
    for(j in i:dimSigma){
	  # compute (d/d j) S^-1
      SigmaiInv.d1j <- Sigmai.inv %*% d1Sigma[,,j]
      dSj <- -SigmaiInv.d1j %*% Sigmai.inv
      dS.j <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSj)
	  # compute (d/di dj) S^-1
      dS.ij <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=d2Sigma[[i]][,,j])
      
	  # compute second derivatives of Sigma^-1 (Harville Ch15, Eq 9.2)
      d2S <- (- Sigmai.inv %*% d2Sigma[[i]][,,j] + 
               SigmaiInv.d1i %*% SigmaiInv.d1j +
               SigmaiInv.d1j %*% SigmaiInv.d1i) %*% Sigmai.inv 
               
      dSij <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=d2S)
      
      marg.hesse[i,j] <- marg.hesse[j,i] <- d2logDet[i,j] -
                        0.5* t(randomEffects) %*% dSij %*% randomEffects -
                        0.5* sum(diag(-F.inv.RE %*% dS.j %*% F.inv.RE %*% dS.i + F.inv.RE %*% dSij))    
    }
    
  }  
  
 marg.Fisher <- as.matrix(-marg.hesse)
   
  return(marg.Fisher)
}


## first and second derivatives of the covariance matrix
dSigma1 <- function(sd,corr){
  derivs <- array(2*exp(2*sd), c(1,1,1)) 
  return(derivs)
}

#d1: result of dSigma1
d2Sigma1 <- function(sd,corr,d1){
  return(list(dsd1=2*d1))
}

dSigma2 <- function(sd,corr){
  derivs <- array(0,c(2,2,3))
  
  dSigma <- diag(2*exp(2*sd))
  if(length(corr)>0){
    dSigma[1,2] <- dSigma[2,1] <- exp(sum(sd[1:2]))*corr[1]/sqrtOf1pr2(corr[1])
    
    # derivative of corr_1
    derivs[2,1,3] <- derivs[1,2,3] <- exp(sum(sd[1:2]))/(sqrtOf1pr2(corr[1])^3)
  }
  
  derivs[,,1:2] <- dSigma
  # derivative of sd_1
  derivs[2,2,1] <- 0
  # derivative of sd_2
  derivs[1,1,2] <- 0

  return(derivs)
}

d2Sigma2 <- function(sd,corr, d1){
  
  derivs <- array(0,c(2,2,3))
  result <- list(dsd1=d1, dsd2=derivs, dcorr1=derivs)
  result$dsd1[1,1,1] <- 2*d1[1,1,1]
  result$dsd1[2,2,2] <- 0
  
  result$dsd2[,,2:3]<- d1[,,2:3]
  result$dsd2[2,2,2] <- 2*d1[2,2,2]
    
    
  if(length(corr)>0){
    result$dcorr1[2,1,3] <- result$dcorr1[1,2,3] <- -(3*corr[1]*exp(sum(sd[1:2])))/(sqrtOf1pr2(corr[1])^5)
  }
  
  return(result)
}


dSigma3 <- function(sd,corr){
  derivs <- array(0,c(3,3,6))
  
  dSigma <- diag(2*exp(2*sd)) #
  if(length(corr)>0){
    dSigma[1,2] <- dSigma[2,1] <- exp(sum(sd[1:2]))*corr[1]/sqrtOf1pr2(corr[1]) #
    dSigma[1,3] <- dSigma[3,1] <- exp(sum(sd[c(1,3)]))*corr[2]/sqrtOf1pr2(corr[2]) #
    dSigma[2,3] <- dSigma[3,2] <- exp(sum(sd[c(2,3)]))*(corr[1]*corr[2]*sqrtOf1pr2(corr[3])+corr[3])/prod(sqrtOf1pr2(corr[1:3]))#
    
    # derivative of corr_1
    derivs[2,1,4] <- derivs[1,2,4] <- exp(sum(sd[1:2]))/(sqrtOf1pr2(corr[1])^3)
    derivs[3,2,4] <- derivs[2,3,4] <-(exp(sum(sd[2:3]))*(corr[2]*sqrtOf1pr2(corr[3])-prod(corr[c(1,3)])))/ (prod(sqrtOf1pr2(corr[2:3]))*(sqrtOf1pr2(corr[1])^3))#
    
    # derivative of corr_2
    derivs[3,1,5] <- derivs[1,3,5] <- exp(sum(sd[c(3,1)]))/(sqrtOf1pr2(corr[2])^3)#
    derivs[3,2,5] <- derivs[2,3,5] <- (exp(sum(sd[2:3]))*(corr[1]*sqrtOf1pr2(corr[3])-prod(corr[c(2,3)])))/ (prod(sqrtOf1pr2(corr[c(1,3)]))*(sqrtOf1pr2(corr[2])^3)) #
    
    # derivative of corr_3
    derivs[3,2,6] <- derivs[2,3,6] <- exp(sum(sd[2:3]))/ (prod(sqrtOf1pr2(corr[c(1,2)]))*(sqrtOf1pr2(corr[3])^3))
  }
  
  derivs[,,1:3] <- dSigma
  # derivative of sd_1
  derivs[2:3,2:3,1] <- 0
  # derivative of sd_2
  derivs[1,c(1,3),2] <- derivs[3,c(1,3),2] <- 0
  # derivative of sd_3
  derivs[1:2,1:2,3] <- 0
  
  return(derivs)

}

d2Sigma3 <- function(sd,corr, d1){
  
  derivs <- array(0,c(3,3,6))
  result <- list(dsd1=d1, dsd2=derivs, dsd3=derivs, dcorr1=derivs, dcorr2=derivs, dcorr3=derivs)
  
  result$dsd1[1,1,1] <- 2*d1[1,1,1]
  result$dsd1[2,2:3,2] <- result$dsd1[3,2,2] <- 0
  result$dsd1[2:3,2:3,3] <-  0 #
  
  result$dsd2[,,2]<- d1[,,2]
  result$dsd2[2,2,2] <- 2*d1[2,2,2]
  result$dsd2[3,2,3] <- result$dsd2[2,3,3] <- d1[3,2,3]#
    
  result$dsd3[,,3]<- d1[,,3]
  result$dsd3[3,3,3] <- 2*d1[3,3,3]#
  
  dSigma <- diag(4*exp(2*sd))
  if(length(corr)>0){
    result$dsd1[2:3,2:3,4] <-  0
    result$dsd1[2:3,2:3,5] <-  0
    result$dsd1[,,6] <-  0
    
    result$dsd2[,,c(4,6)] <- d1[,,c(4,6)]
    result$dsd2[3,2,5] <- result$dsd2[2,3,5] <- d1[3,2,5]
    
    result$dsd3[3,2,4] <- result$dsd3[2,3,4] <- d1[3,2,4]
    result$dsd3[,,c(5,6)] <- d1[,,c(5,6)]
  
    # derivative of corr_1
    result$dcorr1[2,1,4] <- result$dcorr1[1,2,4] <- -(exp(sum(sd[1:2]))*3*corr[1])/(sqrtOf1pr2(corr[1])^5) #
    result$dcorr1[3,2,4] <- result$dcorr1[2,3,4] <- -(exp(sum(sd[2:3]))*(corr[1]*(3*corr[2]*sqrtOf1pr2(corr[3])-2*prod(corr[c(1,3)])) + corr[3]) )/ (prod(sqrtOf1pr2(corr[2:3]))*(sqrtOf1pr2(corr[1])^5)) #
    
    result$dcorr1[3,2,5] <- result$dcorr1[2,3,5] <- (exp(sum(sd[2:3]))*(sqrtOf1pr2(corr[3])+prod(corr[1:3])))/ (prod(sqrtOf1pr2(corr[c(1,2)])^3)*sqrtOf1pr2(corr[3]))

    result$dcorr1[3,2,6] <- result$dcorr1[2,3,6] <- -(exp(sum(sd[2:3]))*corr[1])/ (prod(sqrtOf1pr2(corr[c(1,3)])^3)*sqrtOf1pr2(corr[2]))
      
    # derivative of corr_2
    result$dcorr2[3,1,5] <- result$dcorr2[1,3,5] <- -(exp(sum(sd[c(3,1)]))*3*corr[2])/(sqrtOf1pr2(corr[2])^5)
    result$dcorr2[3,2,5] <- result$dcorr2[2,3,5] <- -(exp(sum(sd[2:3]))*(corr[2]*(3*corr[1]*sqrtOf1pr2(corr[3])-2*prod(corr[c(2,3)])) + corr[3]) )/ (prod(sqrtOf1pr2(corr[c(1,3)]))*(sqrtOf1pr2(corr[2])^5))
      
    result$dcorr2[3,2,6] <- result$dcorr2[2,3,6] <- -(exp(sum(sd[2:3]))*sqrtOf1pr2(corr[2]))/ (prod(sqrtOf1pr2(corr[c(2,3)])^3)*sqrtOf1pr2(corr[1]))
    
    # derivative of corr_3
    result$dcorr3[3,2,6] <- result$dcorr3[2,3,6] <- -(exp(sum(sd[2:3]))*3*corr[3])/ (prod(sqrtOf1pr2(corr[c(1,2)]))*sqrtOf1pr2(corr[3])^5)
    
    
  }
  

  return(result)
}

updateRegression <- function(theta,sd.corr,model=model, 
                             control=list(scoreTol=1e-5, paramTol=1e-8,F.inc=0.01, stepFrac=0.5,niter=30,scale.par=1), 
                             verbose=0, method="nlminb",...){
  regressionParams <- function(theta,...){
    ll <- penLogLik(theta,...)
    attr(ll,"score") <- penScore(theta,...)
    attr(ll,"fisher") <- penFisher(theta,...)
    return(ll)
  }
  regressionParamsMin <- function(theta,...){
    ll <- - penLogLik(theta,...)
    attr(ll,"gradient") <- - penScore(theta,...)
    attr(ll,"hessian") <- penFisher(theta,...)
    return(ll)  
  }
  
  if(any(is.na(sd.corr))){
    cat("WARNING: at least one variance estimates not a number\n")
  } 
  
  if(any(theta == -20)){
    cat("at least one theta reached lower bound -20\n")
  }

  
  lowerBound <- rep(-Inf,model$nFE + model$nOverdisp + model$nRE)
  indexAR <- c(grep("ar.ri",model$namesFE), grep("ar.1",model$namesFE),
               grep("ne.ri",model$namesFE), grep("ne.1",model$namesFE)) 
  lowerBound[indexAR] <- -20
  
  scale <- control$scale.par
  if(is.null(scale)) scale <- 1
  
  # estimate regression coefficients (fixed + random)
  if(method=="nr"){
    res <- newtonRaphson(x=theta,fn=regressionParams,sd.corr=sd.corr,model=model,control=control,verbose=verbose)
    theta.new <- res$coefficients
    ll <- res$loglik
  } else if(method == "nlminb"){
    ll <- function(theta,...) {- penLogLik(theta,...)}
    gr <- function(theta,...) {- penScore(theta,...)}
    trace <- switch(verbose+1,0,0,10,1) # 0-1: no info, 2 every 10th iteration, >=3, each iteration
    if(is.null(trace)) trace <- 1
    res <- nlminb(start=theta, objective=ll, gradient=gr, hessian=penFisher, 
                  sd.corr=sd.corr,model=model, 
                  control=list(trace=trace,abs.tol=1e-20,rel.tol=1e-10,x.tol=1.5e-8, iter.max=control$niter), 
                  lower=lowerBound,
                  scale=scale,...)
    theta.new <- res$par
    ll <- - res$objective
  } else if(method == "nlm"){
    res <- nlm(p=theta, f=regressionParamsMin, sd.corr=sd.corr,model=model,fscale=1,hessian=TRUE,print.level=verbose)
    theta.new <- res$estimate
    ll <- - res$minimum
  } else {
    res <- optim(theta,penLogLik, penScore, sd.corr=sd.corr,model=model,method=method, hessian=TRUE, control=list(fnscale=-1,maxit=1000,trace=verbose))
    theta.new <- res$par
    ll <- res$value
  }
  fisher.unpen <- attributes(penFisher(theta.new,sd.corr=sd.corr,model=model,attributes=TRUE))$fisher
  ll.unpen <- sum(attributes(penLogLik(theta.new,sd.corr=sd.corr,model=model))$loglik)
  
  rel.tol <-  max(abs(theta.new - theta))/max(abs(theta))
  
  if(verbose>0){
    cat("Update for regression parameters:  max|x_0 - x_1| / |x_0|= ", rel.tol,"",ifelse(res$convergence==0,"\n\n", "\nWARNING: algorithm did NOT converge\n\n"))
    print(res$message)
  }  
        
  return(list(par=theta.new, ll=ll, fisher.unpen=fisher.unpen,ll.unpen=ll.unpen, rel.tol=rel.tol, convergence=res$convergence))
}

updateVariance <- function(sd.corr,theta,model, control=list(scoreTol=1e-5, paramTol=1e-8,F.inc=0.01, stepFrac=0.5,niter=30,scale.var=1), verbose=0, method="nlminb", ...){
  # only fixed effects => no variance 
  if(model$nSigma ==0){
    if(verbose>0) cat("No update for variance components\n\n")
    return(list(par=sd.corr,ll=NA,rel.tol=0,convergence=0))
  }
  
  varianceParamsMin <-  function(sd.corr,...){
    ll <- - marLogLik(sd.corr,...)
    attr(ll,"gradient") <- - marScore(sd.corr,...)
    attr(ll,"hessian") <- marFisher(sd.corr,...)
    return(ll)  
  }
  
  
  scale <- control$scale.var
  if(is.null(scale)) scale <- 1


  # estimate variance components
  if(method == "nr"){
    res <- newtonRaphson(x=sd.corr,fn=marLogLik,theta=theta,model=model,control=control, verbose=verbose,...)
    sd.corr.new <- res$coefficients
    ll <- res$loglik
  } else if(method == "nlminb"){
    ll <- function(sd.corr,...) { -c(marLogLik(sd.corr,...))}
    gr <- function(sd.corr,...) { - marScore(sd.corr,...)}
    
    trace <- switch(verbose+1,0,0,10,1) # 0-1: no info, 2 every 10th iteration, >=3, each iteration
    if(is.null(trace)) trace <- 1
    res <- nlminb(start=sd.corr, objective=ll, gradient=gr, hessian=marFisher,
                 theta=theta, model=model,
                 control=list(trace=trace, iter.max=control$niter),scale=scale,...)
    sd.corr.new <- res$par
    ll <- -res$objective
  } else if(method == "nlm"){
    res <- nlm(p=sd.corr, f=varianceParamsMin, theta=theta, model=model,fscale=1,hessian=TRUE, print.level=verbose,...)
    sd.corr.new <- res$estimate
    ll <- -res$minimum
  } else {
    res <- optim(sd.corr,marLogLik, marScore, theta=theta, model=model, method=method, hessian=TRUE, control=list(fnscale=-1,maxit=1000,trace=verbose))
    sd.corr.new <- res$par  
    ll <- res$value

  }
  
  
  if(any(is.na(sd.corr.new))){
    cat("WARNING: at least one variance estimates not a number, no update of variance\n")
    sd.corr.new[is.na(sd.corr.new)] <- -.5
    ll <- c(marLogLik(sd.corr,theta,model))
    return(list(par=sd.corr,ll=ll, rel.tol=NA, convergence=99) ) 
  
  } 
  
  rel.tol <-  max(abs(sd.corr.new - sd.corr))/max(abs(sd.corr))
  
  if(verbose>0)
  cat("Update for variance parameters:  max|x_0 - x_1| / |x_0|= ", rel.tol,"",ifelse(res$convergence==0,"\n\n", "\nWARNING: algorithm did NOT converge\n\n"))
  
  return(list(par=sd.corr.new,ll=ll, rel.tol=rel.tol, convergence=res$convergence) ) 
}



fitHHH <- function(theta,sd.corr,model, control=list(tol=1e-5,niter=15), 
                 cntrl.update=list(scoreTol=1e-5, paramTol=1e-7,F.inc=0.01, stepFrac=0.5,niter=15,
                 scale.par=1, scale.var=1),verbose=0, shrinkage=FALSE, method="nlminb"){
  
  convergence <- 99
  i <- 0
  
  dimFE <- model$nFE + model$nOverdisp
  dimRE <- model$nRE
  
  
  while(convergence != 0 & (i< control$niter)){
    i <- i+1
    # update regression coefficients
    parReg <- updateRegression(theta=theta,sd.corr=sd.corr,model=model, control=cntrl.update,verbose=verbose,method=method)
    
    if(parReg$convergence!=0 & verbose>0)
      cat("Update of regression coefficients in iteration ",i," failed\n")
    if(parReg$convergence >20 && shrinkage){
      theta <- parReg$par
     cat("\n\n***************************************\nshrinkage", 0.1*theta[abs(theta)>10],"\n")
      theta[abs(theta)>10] <- 0.1*theta[abs(theta)>10]
      diag(parReg$fisher.unpen) <- diag(parReg$fisher.unpen)+1e-2
    } else theta <- parReg$par
    
    # update variance components
    parVar <- updateVariance(sd.corr=sd.corr,theta=theta, model=model, fisher.unpen=parReg$fisher.unpen,     
                               control=cntrl.update,verbose=verbose, method=method)
    if(parVar$convergence!=0 & verbose>0)
      cat("Update of variance components in iteration ",i," failed\n")
    #if(parVar$convergence==0)
      sd.corr <- parVar$par
    
    # convergence ?
    if( (parReg$rel.tol < control$tol) && (parVar$rel.tol < control$tol) && (parReg$convergence==0) && (parVar$convergence==0))
      convergence <- 0 
  }
  
  ll <- c(penLogLik(theta=theta,sd.corr=sd.corr,model=model))
  fisher <- penFisher(theta=theta,sd.corr=sd.corr,model=model)
  fisher.var <- marFisher(sd.corr=sd.corr,theta=theta,model=model)
  
  return(list(fixef=head(theta,dimFE),ranef=tail(theta,dimRE),sd.corr=sd.corr, loglik=ll, fisher=fisher, fisherVar=fisher.var, convergence=convergence, dim=c(fixed=dimFE,random=dimRE)))
}


#####################
# x - initial parameter values
# scoreTol - convergence if max(abs(score)) < scoreTol
# paramTol - convergence if rel change in theta < paramTol
# F.inc - eigenvalues of the hessian are computed when the Cholesky factorization
#         fails, and a constant added to the diagonal to make the smallest
#         eigenvalue= F.inc * largest

# fn must return loglikelihood with score and fisher as attributes
#   fn <- function(theta,...){
#     ll <- loglik(theta,...)
#     attr(ll,"score") <- score(theta,...)
#     attr(ll,"fisher") <- fisher(theta,...)
#     return(ll)
#   }
newtonRaphson <- function(x,fn,..., control=list(scoreTol=1e-5, paramTol=1e-8,F.inc=0.01, stepFrac=0.5,niter=30), verbose=FALSE){

  # set default values
  if(is.null(control$scoreTol)) control$scoreTol <- 1e-5
  if(is.null(control$paramTol)) control$paramTol <- 1e-8
  if(is.null(control$F.inc)) control$F.inc <- 0.01
  if(is.null(control$stepFrac)) control$stepFrac <- 0.5
  if(is.null(control$niter)) control$niter <- 30
  
  # number of step reductions, not positive definite Fisher matrices during iterations
  steph <- notpd <- 0
  convergence <- 99
  i <- 0
  x.new <- NA
  
  rel.tol <- function(x,xnew){
    sqrt(sum((xnew-x)^2)/sum(x^2))
  }
  
  score <- function(fn){
    return(attr(fn,"score"))
  }
  fisher <- function(fn){
    return(attr(fn,"fisher"))
  }
  
  ll0 <- c(fn(x,...))
  if(verbose>1) cat("initial loglikelihood",ll0,"\n\n")
  
  # fn cannot be computed at initial par
  if(!is.finite(ll0) | is.na(ll0)){
    cat("fn can not be computed at initial parameter values.\n")
    return(list(convergence=30, notpd = notpd, steph = steph))
  }
      
  while(convergence != 0 & (i< control$niter)){
    i <- i+1
    ll <- fn(x,...)
    if(max(abs(score(ll))) < control$scoreTol){
      convergence <- 0
      break
    }
    # get cholesky decompositon
    F <- fisher(ll)
    F.chol <- try(chol(F),silent=TRUE)
    # could still give a nearly singular matrix
    # => could also check condition number
    
    if(inherits(F.chol,"try-error")){
      if(verbose>1) cat("fisher is not pd\n")
      # fisher is not pd
      notpd <- notpd +1
      ev <- eigen(F,symmetric=TRUE, only.values=TRUE)$values
      #  add a constant to diag(F)
      diag(F) <- diag(F) + (control$F.inc*(max(abs(ev))) - min(ev))/(1-control$F.inc)
      # compute cholesky decomposition of modified fisher
      F.chol <- chol(F)
    }
    direction <- chol2inv(F.chol)%*% score(ll)
    if(max(abs(direction)) < control$paramTol*(max(abs(x))+1e-8) ){
      convergence <- 0
      break
    }
    # do Newton-Raphson step
    x.new <- c(x + direction)
    ll.new <- fn(x.new,...)
    if(verbose>1) cat("iteration",i,"\trel.tol =",rel.tol(x,x.new),"\tabs.tol(score) =",max(abs(score(ll.new))),"\n")
    if(verbose>2) cat("theta =",round(x.new,2),"\n")
    if(verbose>1) cat("loglikelihood =",ll.new,"\n")
      

    ## Backtracking: reduce stepsize until we really improve the loglikelihood
    # ll(x1,lambda) = ll(x0) + lambda * fisher(x0)^-1 %*% score(x0)
    i.backstep <- 0
    ## Gray (2001) Ch 3: Unconstrained Optimization and Solving Nonlinear Equations
    # It is technically possible to construct sequences where ll(x1) > ll(x0)
    # at each step but where the sequence never converges.
    # For this reason a slightly stronger condition is usually used.
    # Dennis and Schnabel (1983): Numerical Methods for Unconstrained 
    #     Optimization and Nonlinear Equations. SIAM. (ch 6,3.2, p.126)
    # recommend requiring that lambda satisfy
    #   ll(x1) > ll(x0) + 1e-4 *(x1-x0)' %*% score(x0)
    while((is.na(ll.new) || (ll.new < c(ll)+ (1e-4)*sum(direction*score(ll)))) & (i.backstep <= 20)){
      if(verbose>1 & i.backstep==0) cat("backtracking: ")
      i.backstep <- i.backstep +1
      steph <- steph +1
      # reduce stepsize by a fixed fraction stepFrac
      direction <- control$stepFrac*direction
      x.new <- c(x + direction)
      ll.new <- fn(x.new,...)
      if(verbose>1) cat("*")
    }
    if(verbose & i.backstep>0) cat("\n")
    if(i.backstep >20){
      if(verbose>1)cat("backtracking did not improve fn\n")
      #cat("ll: ",ll,"\tll.new: ",ll.new,"\n")
      convergence <- 10
      break
    }

    x <- c(x.new)
      
    if(verbose>1) cat("\n")
  }
  ll <- fn(x,...)
  
  # max number of iterations reached, but check for convergence
  if(max(abs(score(ll))) < control$scoreTol){
    convergence <- 0
  }

  # convergence if
  # 1) relative difference between parameters is small
  # 2) absolute value of gradient is small
  # 3) stop after niter iterations
  
  if(i==control$niter & convergence !=0){
      if(verbose>1) cat("Newton-Raphson stopped after",i,"iterations!\n")
      # iteration limit reached without convergence
      convergence <- 10
  }
  if(verbose>1) cat("iteration",i,"\trel.tol =",rel.tol(x,x.new),"\tabs.tol(score) =",max(abs(score(ll))),"\n")
  if(verbose>2) cat("theta =",round(x.new,2),"\n")
  if(verbose>1) cat("loglikelihood =",c(ll),"\n\n")
 
  # loglikelihood
  loglik <- c(ll)
  
  # fisher info
  F <- fisher(ll)
  if(inherits(try(solve(F),silent=TRUE),"try-error")){ 
    cat("\n\n***************************************\nfisher not regular!\n")
    #print(summary(x))
    return(list(coefficients=x, loglikelihood=loglik, fisher=FALSE, convergence=22, notpd = notpd, steph = steph))
  }
  
  # check if solution is a maximum (i.e. if fisher is pd )
  eps <- 1e-10
  if(!all(eigen(F,symmetric=TRUE, only.values=TRUE)$values > eps)){
    if(verbose>1) cat("fisher information at solution is not pd\n")
    return(list(coefficients=x, loglikelihood=loglik, fisher=FALSE, convergence=21, notpd = notpd, steph = steph))
  }
 
  if(verbose>0) 
    cat("number of iterations = ",i," coverged = ", convergence ==0," log-likelihood = ",loglik, " notpd = ", notpd, " steph = ", steph, "\n")
  result <- list(coefficients=x, loglikelihood=loglik, fisher=FALSE, 
                 convergence=convergence, notpd=notpd, steph=steph,niter=i)
  return(result)
  
}

##############
addSeason2formula <- function(f=~1,       # formula to start with
                              S=1,         # number of sine/cosine pairs
                              period=52
                              ){
  # return formula as is if S = 0
  if(max(S) == 0) return(f)
  
  f <- deparse(f)
  # create formula
  if(max(S)>0 & length(S)==1){
    for(i in 1:S){
      f <- paste(f,"+sin(",2*i,"*pi*t/",period,")+cos(",2*i,"*pi*t/",period,")",sep="")
    }
  } else {
    nSeason <- length(S)
    for(i in 1:max(S)){
      which <- paste(i <= S,collapse=",")
      f <- paste(f,"+ fe( sin(",2*i,"*pi*t/",period,"), which=c(",which,")) + fe( cos(",2*i,"*pi*t/",period,"), which=c(",which,"))",sep="")
    }
  }
  return(as.formula(f))
}


oneStepAhead <- function(result, # result of call to hhh4
                         tp,     # one-step-ahead predictions for time points (tp+1):nrow(stsObj)
                         verbose=TRUE,
                         keep.estimates = FALSE){
  
  stsObj <- result$stsObj
  control <- result$control
  if(is.null(result$terms)){
	# get model terms
	model <- interpretControl(control,stsObj)
  } else {
	model <- result$terms
  }
  
  # which model: negbin or poisson?
  dimOverdisp <- model$nOverdisp
  negbin <- dimOverdisp>0

  nTime <- nrow(stsObj)
  pred <- matrix(NA,nrow=length(tp:nTime)-1,ncol=ncol(stsObj))
  psi <- matrix(NA,nrow=length(tp:nTime)-1,ncol=ifelse(dimOverdisp>1,ncol(stsObj),1))
    
  res <- resN <- result
  coefs <- coef(res,reparamPsi=FALSE)
  control.i<-control
  
  which <- 1
  search2 <- function(i,which){
    control.i$subset <- 2:i
    
    if(which==1){
      #starting values of previous time point i-1
      control.i$start <- list(fixed=fixef(res.old),random=ranef(res.old),sd.corr=getSdCorr(res.old))
    } else {
      #starting values of last time point
      control.i$start <- list(fixed=fixef(resN),random=ranef(resN),sd.corr=getSdCorr(resN))
    }
    res <- hhh4(stsObj,control.i)  

    return(res)
  }
  
  if(keep.estimates){
	# save value of log-likelihood (pen+mar) and parameter estimates
	params <- matrix(NA,nrow=length(tp:nTime)-1,ncol=length(coef(res))+2)
	# save values of estimated covariance
	vars <- matrix(NA,nrow=length(tp:nTime)-1,ncol=model$nSigma)
  } else {
	params <- vars <- NULL
  }

  # make one-step ahead prediction for time point tp+1
  for(i in (tp):(nTime-1)) {
      cat(nTime-i,"\n")
      res.old <- res
      
      # fit model to data for time points 1,...,tp
      # use initial values from previous fit or from fit to all
      res <- search2(i, which=which)
      
      # check convergence        
      # do a grid search in case of non-convergence ?
      if(res$convergence){
        # make one-step ahead prediction for time point tp+1
        pred[i-tp+1,] <- tail(predict(res,newSubset=2:(i+1),type="response"),1)
		# get overdispersion parameter
        if(negbin)
          psi[i-tp+1,] <- coef(res, reparamPsi=FALSE)[grep("overdisp",names(coef(res)))]
		# save parameter estimates
		if(keep.estimates){
		  params[i-tp+1,] <- c(res$loglikelihood,res$margll,coef(res,reparamPsi=FALSE))
		  vars[i-tp+1,] <- getSdCorr(res)
		}
      } else {
        cat("NO convergence in iteration", i-tp,"\n")
        res <- res.old
      }
  }
   
  return(list(mean=pred, psi=psi, x=observed(stsObj[(tp+1):nTime]),
		      allConverged=all(!is.na(pred)),
              params=params,variances=vars))

}

