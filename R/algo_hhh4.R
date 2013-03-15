################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### hhh4 is an extended version of algo.hhh for the sts-class
### The function allows the incorporation of random effects and covariates.
###
### Copyright (C) 2010-2013 Michaela Paul and Sebastian Meyer
### $Revision: 483 $
### $Date: 2013-01-24 16:36:51 +0100 (Do, 24. Jan 2013) $
################################################################################

# - some function arguments are currently not used (but will eventually)
# - formula formulation is experimental and not yet implemented in full generality 
# - do some profiling...

## Error message issued in loglik, score and fisher functions upon NA parameters
ADVICEONERROR <- "\n  Try different starting values or another optimizer.\n"


### Default control argument for hhh4
### This argument list is also set as the formal 'control' argument of hhh4()
### below the definition of hhh4

CONTROL.hhh4 <- alist(
    ar = list(f = ~ -1,        # a formula " exp(x'lamba)*y_t-1 "
                               # (ToDo: or a matrix " Lambda %*% y_t-1 ")
              lag = 1,         # autoregression on y_i,t-lag (currently unused)
              weights = NULL,  # for a contact matrix model (currently unused)
              initial = NULL   # initial parameters values if pred is a matrix
                               # or if pred = ~1 (not really used ATM)
    ),
    ne = list(f = ~ -1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
              lag = 1,         # regression on y_j,t-lag  (currently unused)
              weights = neighbourhood(stsObj),  # weights w_ji
              initial = NULL   # initial parameter values if pred = ~1 (not really used ATM)
    ),
    end = list(f = ~ 1,        # a formula " exp(x'nu) * n_it "
               offset = 1,     # optional offset e_it
               initial = NULL  # initial parameter values if pred = ~1 (not really used ATM)
    ),
    family = c("Poisson", "NegBin1", "NegBinM"),
    subset = 2:nrow(stsObj),   # typically 2:nTime if model contains lags
    optimizer = list(stop = list(tol=1e-5, niter=100),   # control arguments
                     regression = list(method="nlminb"), # for optimization 
                     variance = list(method="nlminb")),
    verbose = FALSE,           # level of reporting during hhh4() processing
    start = list(fixed=NULL,random=NULL,sd.corr=NULL),  # list with initials,
                               # overrides any initial values in formulas 
    data = data.frame(t=epoch(stsObj)-1), # data.frame or named list with
                               # covariates that appear in any of the formulae
    keep.terms = FALSE
)



### Main function, which is to be called by the user

hhh4 <- function (stsObj, control, check.analyticals = FALSE)
{
  ## Convert old disProg class to new sts class
  if(inherits(stsObj, "disProg")) stsObj <- disProg2sts(stsObj)
  
  ## check control and set default values (for missing arguments)
  control <- setControl(control, stsObj)
  #* if univariate, no neighbouring term
  #* check nhood matrix
  
  ## get model terms
  model <- interpretControl(control, stsObj)
  dimFixedEffects <- model$nFE + model$nd + model$nOverdisp
  dimRandomEffects <- model$nRE

  ## starting values 
  #* -> better default values possible
  theta.start <- model$initialTheta
  Sigma.start <- model$initialSigma
  
  ## check if initial values are valid
  ## CAVE: there might be NA's in mu if there are missing values in Y
  mu <- meanHHH(theta.start,model)$mean
  if(any(mu==0, na.rm=TRUE) || any(is.infinite(mu)))
    stop("some mean is degenerate (0 or Inf) at initial values")

  ## check score vector and fisher information at starting values
  check.analyticals <- if (isTRUE(check.analyticals)) {
      if (length(theta.start) > 50) "maxLik" else "numDeriv"
  } else if (is.character(check.analyticals)) {
      match.arg(check.analyticals, c("numDeriv", "maxLik"), several.ok=TRUE)
  } else NULL
  if (length(check.analyticals) > 0L) {
      cat("\nPenalized log-likelihood:\n")
      resCheckPen <- sapply(check.analyticals, function(derivMethod) {
          if (require(derivMethod, character.only=TRUE)) {
              do.call(paste("checkDerivatives", derivMethod, sep="."),
                      args=alist(penLogLik, penScore, penFisher, theta.start,
                      sd.corr=Sigma.start, model=model))
          }
      }, simplify=FALSE, USE.NAMES=TRUE)
      if (length(resCheckPen) == 1L) resCheckPen <- resCheckPen[[1L]]
      resCheckMar <- if (length(Sigma.start) == 0L) list() else {
          cat("\nMarginal log-likelihood:\n")
          fisher.unpen <- attr(penFisher(theta.start, Sigma.start, model,
                                         attributes=TRUE), "fisher")
          resCheckMar <- sapply(check.analyticals, function(derivMethod) {
              if (require(derivMethod, character.only=TRUE)) {
                  do.call(paste("checkDerivatives", derivMethod, sep="."),
                          args=alist(marLogLik, marScore, marFisher, Sigma.start,
                          theta=theta.start, model=model,
                          fisher.unpen=fisher.unpen))
              }
          }, simplify=FALSE, USE.NAMES=TRUE)
          if (length(resCheckMar) == 1L) resCheckMar[[1L]] else resCheckMar
      }
      resCheck <- list(pen = resCheckPen, mar = resCheckMar)
      return(resCheck)
  }

  ## maximize loglikelihood (penalized and marginal)
  myoptim <- fitHHH(theta=theta.start,sd.corr=Sigma.start, model=model,
                    cntrl.stop       = control$optimizer$stop,
                    cntrl.regression = control$optimizer$regression,
                    cntrl.variance   = control$optimizer$variance,
                    verbose=control$verbose)

  convergence <- myoptim$convergence == 0
  thetahat <- c(myoptim$fixef, myoptim$ranef)
  loglik <- myoptim$loglik
  fisher <- myoptim$fisher
  cov <- try(solve(fisher), silent=TRUE) # approximation to inverted fisher info

  ## check for degenerate fisher info
  if(inherits(cov, "try-error")){ # fisher info is singular
    cat("Fisher info singular!\n")
    convergence <- FALSE
  } else if(any(!is.finite(diag(cov))) || any(diag(cov)<0)){
    cat("Infinite or negative cov!\n")
    convergence <- FALSE
  }

  if (!convergence) {
      cat("Penalized loglikelihood =", loglik, "\n")
      thetastring <- paste(round(thetahat,2), collapse=", ")
      thetastring <- strwrap(thetastring, exdent=10, prefix="\n", initial="")
      cat("theta = (", thetastring, ")\n")
      cat("WARNING: Results are not reliable!", ADVICEONERROR)
      res <- myoptim
      res$convergence <- convergence
      res$call <- match.call()
      class(res) <- "ah4"
      return(res)
  }

  ## optimization successful, return a full "ah4" object
  fitted <- meanHHH(thetahat,model)$mean

  if(dimRandomEffects>0){
    Sigma.orig <- myoptim$sd.corr
    Sigma.cov <- solve(myoptim$fisherVar)
    dimnames(Sigma.cov) <- list(names(Sigma.orig),names(Sigma.orig))
  } else {
    Sigma.orig <- Sigma.cov <- NULL
  }

  se <- sqrt(diag(cov))

  margll <- marLogLik(Sigma.orig, thetahat, model)
  Sigma.trans <- getSigmai(head(Sigma.orig,model$nVar),
                           tail(Sigma.orig,model$nCorr),
                           model$nVar)

  ## Done
  result <- list(coefficients=thetahat, se=se, cov=cov, 
                 Sigma=Sigma.trans,     # estimated covariance matrix of ri's
                 Sigma.orig=Sigma.orig, # variance parameters on original scale
                 Sigma.cov=Sigma.cov,   # covariance matrix of Sigma.orig
                 call=match.call(),
                 dim=c(fixed=dimFixedEffects,random=dimRandomEffects),
                 loglikelihood=loglik, margll=margll, 
                 convergence=convergence,
                 fitted.values=fitted,
                 control=control,
                 terms=if(control$keep.terms) model else NULL,
                 stsObj=stsObj,         # FIXME: stsObj also in 'control$data'
                 lag=1, nObs=sum(!model$isNA[control$subset,]),
                 nTime=length(model$subset), nUnit=ncol(stsObj))
  class(result) <- "ah4"
  result
}

## Store default control arguments in the formal argument list of hhh4
formals(hhh4)[["control"]] <- as.pairlist(CONTROL.hhh4)


## set default values for model specifications in control
setControl <- function (control, stsObj)
{
  stopifnot(is.list(control))
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)
  if(nTime <= 2) stop("too few observations")
  
  ## arguments in 'control' override any corresponding default arguments
  defaultControl <- lapply(CONTROL.hhh4, eval, envir=environment())
  environment(defaultControl$ar$f) <- environment(defaultControl$ne$f) <-
      environment(defaultControl$end$f) <- .GlobalEnv
  defaultControl$data <- as.list(defaultControl$data) # since we will add stsObj
  ##<- NROW(stsObj) only works as intended since R 2.15.0, but this is required
  ##   for adding stsObj to a data.frame; thus we switch to a list
  control <- modifyList(defaultControl, control)

  ## add stsObj (maybe overwrite an old one from the control object)
  control$data$.sts <- stsObj
  
  ## add nTime and nUnits to control list
  control$nTime <- nTime                # FIXME: are these actually
  control$nUnits <- nUnits              #        used anywhere?

  ## check that component specifications are list objects
  for (comp in c("ar", "ne", "end")) {
      if(!is.list(control[[comp]])) stop("'control$", comp, "' must be a list")
  }

  ## check lags in "ar" and "ne" components
  for (comp in c("ar", "ne")) {
      if (!isScalar(control[[comp]]$lag))
          stop("'control$", comp, "$lag' must be scalar")
      control[[comp]]$lag <- as.integer(control[[comp]]$lag)
      if (control[[comp]]$lag != 1L)
          stop("currently, only 'control$", comp, "$lag=1' is supported")
  }

  
  ### check AutoRegressive component
  
  control$ar$isMatrix <- is.matrix(control$ar$f)

  if (is.matrix(control$ar$f)) {
      if (any(dim(control$ar$f) != nUnits))
          stop("'control$ar$f' must be a square matrix of size ", nUnits)
      # use identity matrix if weight matrix is missing
      if (is.null(control$ar$weights)) control$ar$weights <- diag(nrow=nUnits)
      control$ar$inModel <- TRUE
  } else if (inherits(control$ar$f, "formula")) {
      if (!is.null(control$ar$weights)) {
          warning("argument 'control$ar$weights' is not used")
          control$ar$weights <- NULL
      }
      # check if formula is valid
      control$ar$inModel <- isInModel(control$ar$f)
  } else {
      stop("'control$ar$f' must be either a formula or a matrix")
  }

  if (is.matrix(control$ar$weights)) {
      if (any(dim(control$ar$weights) != nUnits))
          stop("'control$ar$weights' must be a square matrix of size ", nUnits)
  }
  
  
  ### check NEighbourhood component
  
  if (!inherits(control$ne$f, "formula"))
      stop("'control$ne$f' must be a formula")
  control$ne$inModel <- isInModel(control$ne$f)
  
  if (control$ne$inModel) {
      ## if ar$f is a matrix it includes neighbouring units => no "ne" component
      if (control$ar$isMatrix)
          stop("there must not be an extra \"ne\" component ",
               "if 'control$ar$f' is a matrix")
      ## check ne$weights specification
      checkWeights(control$ne$weights, nUnits, nTime,
                   neighbourhood(stsObj), control$data)
  } else control$ne$weights <- NULL

  
  ### check ENDemic component
  
  if (!inherits(control$end$f, "formula"))
      stop("'control$end$f' must be a formula")
  control$end$inModel <- isInModel(control$end$f)

  if (is.matrix(control$end$offset) && is.numeric(control$end$offset)){
      if (!identical(dim(control$end$offset), dim(stsObj)))
          stop("'control$end$offset' must be a numeric matrix of size ",
               nTime, "x", nUnits)
      if (any(is.na(control$end$offset)))
          stop("'control$end$offset' must not contain NA values")
  } else if (!identical(as.numeric(control$end$offset), 1)) {
      stop("'control$end$offset' must either be 1 or a numeric ",
           nTime, "x", nUnits, " matrix")
  }


  ### check remaining components of the control list
  
  control$family <- match.arg(control$family, defaultControl$family)  

  if (!is.vector(control$subset, mode="numeric") ||
      !all(control$subset %in% seq_len(nTime)))
      stop("'control$subset' must be %in% 1:", nTime)

  if (!is.list(control$optimizer) ||
      any(! sapply(c("stop", "regression", "variance"),
                   function(x) is.list(control$optimizer[[x]]))))
      stop("'control$optimizer' must be a list of lists")

  if (!length(control$verbose) == 1L ||
      (!is.logical(control$verbose) && !is.numeric(control$verbose)))
      stop("'control$verbose' must be a numeric or logical value")
  control$verbose <- as.integer(control$verbose)

  if (!is.list(control$start) ||
      any(! sapply(c("fixed", "random", "sd.corr"),
                   function(x) is.null(control$start[[x]]) ||
                               is.numeric(control$start[[x]]))))
      stop("'control$start' must be a list of numeric start values")
  
  stopifnot(length(control$keep.terms) == 1L, is.logical(control$keep.terms))

  ## Done
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
    initial <- rep.int(0,dim.fe)
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
    
  } else if(type=="car"){ # construct penalty matrix K
    K <- neighbourhood(sts)
    checkNeighbourhood(K)
    K <- K == 1                         # indicate first-order neighbours
    ne <- colSums(K)                    # number of first-order neighbours
    K <- -1*K
    diag(K) <- ne
    
    dimK <- nrow(K)
    
    # check rank of the nhood, only connected neighbourhoods are allowed
    if(qr(K)$rank != dimK-1) stop("neighbourhood matrix contains islands")
    # singular-value decomposition of K
    svdK <- svd(K)
    # just use the positive eigenvalues  of K in descending order 
    # for a the factorisation of the penalty matrix K = LL'
    L <- svdK$u[,-dimK] %*% diag(sqrt(svdK$d[-dimK]))            #* only use non-zero eigenvalues
    
    # Z = L(L'L)^-1, which can't be simplified to Z=(L')^-1 as L is not square
    Z <- L %*% solve(t(L)%*%L)
    
    dim.re <- dimK-1
    mult <- "%*%"
  } 
  
  #* check length of initial values + set default values?
  if(!is.null(initial.re) & length(initial.re) !=dim.re){
    stop("Initial values for \'",type,"\' must be of length ",dim.re)
  } else if(is.null(initial.re)){
    initial.re <- rnorm(dim.re,0,sd=sqrt(0.001))
  }
  if(!is.null(initial.var) & length(initial.var) !=1){
    stop("Initial values for \'",type,"\' must be of length ",1)
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
  # (only if there are variables in addition to an intercept "1")
  if(nVars >0){
	fe.raw <- setdiff(1:nVars, unlist(attr(term, "specials")))
  } else {
	fe.raw <- numeric(0)
  }
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
  
  res <- cbind(res, deparse.level=0) # ensure res has matrix dimensions
  if(sum(unlist(res["intercept",])) > 1) {
    stop("There can only be one intercept in the formula ", deparse(substitute(f)))
  }
  
  # random effects
  RI <- attr(term, "specials")$ri

  if(length(RI)>1){
    stop("There can only be one random intercept in the formula ",deparse(substitute(f)))
  }
  if(!is.null(RI)){
    for(i in RI){
      res <- cbind(res,c(eval(vars[[i+1]], envir=env),list(offsetComp=component)))
    }      
  }
  
  return(res)
}


## Create function (pars, type = "response") which
## returns the weighted sum of time-lagged counts of neighbours
## (or its derivates, if type = "gradient" or type = "hessian").
## For type="reponse", the result is a matrix/array, otherwise a list of such
## objects, which for the gradient has length length(pars) and
## length(pars)*(length(pars)+1)/2 for the hessian.
neOffsetFUN <- function (Y, neweights, nbmat, data, lag)
{
    if (is.null(neweights)) { # no neighbourhood component
        return(as.function(alist(...=, 0), envir=.GlobalEnv))
    }

    ## return function (living in THIS environment)
    if (is.list(neweights)) { # parametric weights
        function (pars, type = "response") {
            name <- switch(type, response="w", gradient="dw", hessian="d2w")
            weights <- neweights[[name]](pars, nbmat, data)
            ## gradient and hessian are lists if length(pars$d) > 1L
            ## and single matrices/arrays if == 1 => _c_onditional lapply
            res <- clapply(weights, function(W) weightedSumNE(Y, W, lag))
            ##<- clapply always returns a list (possibly of length 1)
            if (type=="response") res[[1L]] else res
        }
    } else { # fixed (known) array (0-length pars)
        initoffset <- weightedSumNE(Y, neweights, lag)
        function (pars, type = "response") initoffset
        ## this will not be called for other types
    }
}


# interpret and check the specifications of each component
# control must contain all arguments, i.e. setControl was used
interpretControl <- function(control, stsObj){

  # get environment for evaluation of covariates
  env <- control$data
  
  nTime <- nrow(stsObj)
  nUnits <- ncol(stsObj)

  Y <- observed(stsObj)

  
  ##########################################################################
  ##  get the model specifications for each of the three components
  ##########################################################################
  ar <- control$ar
  ne <- control$ne
  end <- control$end

  ## FIXME: really weird... if the array below is evaluated, code will be faster
  ## A complex example model takes 132s (without this nonsense array evaluation)
  ## or 116s (with) for the whole fit
  if (ne$inModel && is.array(ne$weights)) {
    dev.null <- array(ne$weights, c(nUnits,nUnits,nTime))
    ## for instance, str(ne$weights), as.vector(ne$weights), or
    ## array(ne$weights, dim(ne$weights)) do not have this effect...
  }

  ## create list of offsets of the three components
  Ym1 <- rbind(matrix(NA_integer_, ar$lag, nUnits), head(Y, nTime-ar$lag))
  Ym1.ne <- neOffsetFUN(Y, ne$weights, neighbourhood(stsObj), control$data,
                        ne$lag)
  offsets <- list(ar=Ym1, ne=Ym1.ne, end=end$offset)

  ## Ym1.ne is a function (pars, type="response)!
  ## Initial parameter vector 'd' of the neighbourhood weight function
  initial.d <- if (is.list(ne$weights)) ne$weights$initial else numeric(0L)
  dim.d <- length(initial.d)
  names.d <- if (dim.d == 0L) character(0L) else {
      paste0("neweights.", if (is.null(names(initial.d))) {
          if (dim.d==1L) "d" else paste0("d", seq_len(dim.d))
      } else names(initial.d))
  }

  ## determine all NA's (offset[[i]] is either 0 or a nTime x nUnits matrix)
  isNA <- (is.na(Ym1) | is.na(Ym1.ne(initial.d)) | is.na(Y))


  ## get terms for all components
  all.term <- NULL
  if(ar$isMatrix){                      # ar$f is a matrix
      stop("not yet implemented")
  } else if(ar$inModel){                # ar$f is a formula
      all.term <- cbind(all.term, checkFormula(ar$f, env=env, component=1))
  }
  if(ne$inModel){
      all.term <- cbind(all.term, checkFormula(ne$f, env=env, component=2))
  }
  if(end$inModel){
      all.term <- cbind(all.term, checkFormula(end$f,env=env, component=3))
  }
  
  dim.fe <- sum(unlist(all.term["dim.fe",]))
  dim.re.group <- unlist(all.term["dim.re",], use.names=FALSE)
  dim.re <- sum(dim.re.group)
  dim.var <- sum(unlist(all.term["dim.var",]))
  dim.corr <- sum(unlist(all.term["corr",]))
  
  if(dim.corr>0){
    if(dim.var!=dim.corr) stop("Use corr=\'all\' or corr=\'none\' ")
    dim.corr <- switch(dim.corr,0,1,3) 
  }
  
  # the vector with dims of the random effects must be equal if they are correlated
  if(length(unique(dim.re.group[dim.re.group>0]))!=1 & dim.corr>0){ 
    stop("Correlated effects must have same penalty")
  }
  
  n <- c("ar","ne","end")[unlist(all.term["offsetComp",])]
  names.fe <- names.var <- NULL
  for(i in 1:length(n)){
    names.fe <- c(names.fe,paste(n[i],all.term["name",i][[1]], sep="."))
    if(all.term["random",i][[1]]) names.var <- c(names.var,paste("sd",n[i],all.term["name",i][[1]], sep="."))
  }
  index.fe <- rep(1:ncol(all.term), times=unlist(all.term["dim.fe",]))
  index.re <- rep(1:ncol(all.term), times=unlist(all.term["dim.re",])) 
  
  # poisson or negbin model
  fam <- match.arg(control$family, c("Poisson","NegBin1","NegBinM"))
  if(fam=="Poisson"){
    ddistr <- function(y,mu,size){
        dpois(y, lambda=mu, log=TRUE)
    }
    dim.overdisp <- 0L
  } else {
    ddistr <- function(y,mu,size){
        dnbinom(y, mu=mu, size=size, log=TRUE)
    }
    ## version that can handle size = Inf (i.e. the Poisson special case):
    ## ddistr <- function (y,mu,size) {
    ##     poisidx <- is.infinite(size)
    ##     res <- y
    ##     res[poisidx] <- dpois(y[poisidx], lambda=mu[poisidx], log=TRUE)
    ##     res[!poisidx] <- dnbinom(y[!poisidx], mu=mu[!poisidx],
    ##                              size=size[!poisidx], log=TRUE)
    ##     res
    ## }
    dim.overdisp <- if (fam=="NegBin1") 1L else nUnits
  }
  environment(ddistr) <- .GlobalEnv     # function is self-contained

  # parameter start values
  initial.fe <- unlist(all.term["initial.fe",])
  initial.overdisp <- rep.int(2, dim.overdisp)
  initial.fe.d.overdisp <- c(initial.fe, initial.d, initial.overdisp)
  initial.re <- unlist(all.term["initial.re",])
  initial.sd <- unlist(all.term["initial.var",])
  initial.corr <- rep.int(0, dim.corr)
  initial.sd.corr <- c(initial.sd, initial.corr)
  
  # check if a list of 'start' values is supplied
  if(!is.null(control$start$fixed)){
    if(length(control$start$fixed) != dim.fe + dim.d + dim.overdisp)
    stop("initial values in start$fixed must be of length ",
         dim.fe + dim.d + dim.overdisp)
    initial.fe.d.overdisp <- control$start$fixed
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

  # set names of parameter vectors
  names.overdisp <- if(dim.overdisp > 1L){
    paste(paste("-log(overdisp", colnames(stsObj), sep=".") ,")", sep="")
  } else {
    rep("-log(overdisp)",dim.overdisp)  # dim.overdisp may be 0
  }
  names(initial.fe.d.overdisp) <- c(names.fe, names.d, names.overdisp)
  names(initial.sd.corr) <- c(names.var, head(paste("corr",1:3,sep="."),dim.corr))

  # Done
  result <- list(response = Y,
                 terms = all.term,
                 nTime = nTime,
                 nUnits = nUnits,
                 nFE = dim.fe,
                 nd = dim.d,
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
                 initialTheta = c(initial.fe.d.overdisp, initial.re),
                 initialSigma = initial.sd.corr,
                 offset = offsets,
                 family = ddistr,
                 subset = control$subset,
                 isNA = isNA
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
  rownames(corr) <- colnames(corr) <-
      sub("^sd\\.","",names(x$Sigma.orig))[grep("^sd\\.",names(x$Sigma.orig))]
  return(corr)
}


splitParams <- function(theta, model){
  fixed <- theta[seq_len(model$nFE)]
  d <- theta[model$nFE + seq_len(model$nd)]
  overdisp <- theta[model$nFE + model$nd + seq_len(model$nOverdisp)]
  random <- theta[seq.int(to=length(theta), length.out=model$nRE)]
  list(fixed=fixed, random=random, overdisp=overdisp, d=d)
}


### compute predictor

meanHHH <- function(theta, model)
{
  # unpack theta
  pars <- splitParams(theta, model)
  fixed <- pars$fixed
  random <- pars$random
  d <- pars$d

  # unpack model
  term <- model$terms
  offsets <- model$offset  # offsets[[i]] is either 0 or a nTime x nUnits matrix
                           # however, offsets[[2]] is a function of 'd'
  nGroups <- model$nGroups
  comp <- unlist(term["offsetComp",])
  idxFE <- model$indexFE
  idxRE <- model$indexRE
  
  subset <- model$subset
  isNA <- model$isNA                    # set missing values in observed Y to NA

  toMatrix <- function(par, r=model$nTime, c=model$nUnits){
    matrix(par,r,c, byrow=TRUE)
  }
  
  # go through groups of parameters and compute the predictor of each component
  computePartMean <- function (component)
  {
    pred <- nullMatrix <- toMatrix(0)
    is.na(pred) <- isNA
    
    if(!any(comp==component)) return(pred)
    
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
      pred <- pred + X*fe + Z.re
    }
    
    exp(pred)                  # CAVE: this is without the multiplicative offset
  }
  
  ## autoregressive component
  ar.mean <- (computePartMean(1) * offsets[[1L]])[subset,,drop=FALSE]
  
  ## neighbourhood component
  ne.exppred <- computePartMean(2)
  ne.offset <- offsets[[2L]](d)         # this is just 0 if no "ne" in model
  ne.mean <- (ne.exppred * ne.offset)[subset,,drop=FALSE]
  
  ## endemic component
  end.mean <- (computePartMean(3) * offsets[[3L]])[subset,,drop=FALSE]
   
  ## total mean
  mean <- ar.mean + ne.mean + end.mean

  ## Done (FIXME: "epidemic" seems to be unused)
  list(mean=mean, epidemic=ar.mean+ne.mean, endemic=end.mean,
       epi.own=ar.mean, epi.neighbours=ne.mean,
       ne.exppred=ne.exppred[subset,,drop=FALSE])
}


############################################


penLogLik <- function(theta, sd.corr, model, attributes=FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE

  ## unpack parameters
  pars <- splitParams(theta, model)
  psi <- exp(pars$overdisp)             # = 1/psi, pars$overdisp = -log(psi)
  if(dimPsi > 1L) {
      psi <- matrix(psi, ncol=model$nUnits, nrow=model$nTime,
                    byrow=TRUE)[subset,,drop=FALSE]
  }
  if (dimRE > 0) {
      randomEffects <- pars$random
      sd   <- head(sd.corr, model$nVar)
      corr <- tail(sd.corr, model$nCorr)
      dimBlock <- model$rankRE[model$rankRE>0]
      Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }
  
  ############################################################
  
  #psi might be numerically equal to 0 or Inf in which cases dnbinom (in meanHHH)
  #would return NaN (with a warning). The case size=Inf rarely happens and
  #corresponds to a Poisson distribution. Currently this case is not handled
  #in order to have the usual non-degenerate case operate faster.
  #For size=0, log(dnbinom) equals -Inf for positive x or if (x=0 and mu=0), and
  #zero if x=0 and mu>0 and mu<Inf. Thus, if there is at least one case in Y
  #(x>0, which is always true), we have that sum(ll.units) = -Inf, hence: 
  if (any(psi == 0)) return(-Inf)

  ## evaluate mean
  mu <- meanHHH(theta, model)
  meanTotal <- mu$mean
  # if, numerically, meanTotal=Inf, log(dnbinom) or log(dpois) both equal -Inf, hence:
  #if (any(is.infinite(meanTotal))) return(-Inf)
  # however, since meanTotal=Inf does not produce warnings below and this is a rare
  # case, it is faster to not include this conditional expression

  ## penalization term for random effects
  lpen <- if (dimRE==0) 0 else { # there are random effects
    ##-.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
    ## the following implementation takes ~85% less computing time !
    -0.5 * crossprod(randomEffects, Sigma.inv) %*% randomEffects
  }
  lpen <- c(lpen)                       # drop 1x1 matrix dimensions
  
  ## log-likelihood
  ddistr <- model$family
  ll.units <- colSums(ddistr(Y,meanTotal,psi), na.rm=TRUE)

  ## penalized log-likelihood
  ll <- sum(ll.units) + lpen

  ## Done
  if (attributes) {
      attr(ll, "loglik") <- ll.units
      attr(ll, "logpen") <- lpen
  }
  ll
}


penScore <- function(theta, sd.corr, model)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE

  ## unpack parameters
  pars <- splitParams(theta, model)
  psi <- exp(pars$overdisp)             # = 1/psi, pars$overdisp = -log(psi)
  if(dimPsi > 1L) {
      psi <- matrix(psi, ncol=model$nUnits, nrow=model$nTime,
                    byrow=TRUE)[subset,,drop=FALSE]
  }
  if (dimRE > 0) {
      randomEffects <- pars$random
      sd   <- head(sd.corr, model$nVar)
      corr <- tail(sd.corr, model$nCorr)
      dimBlock <- model$rankRE[model$rankRE>0]
      Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ## further model unpacking
  term <- model$terms
  nGroups <- model$nGroups
  dimd <- model$nd

  ## evaluate mean
  mu <- meanHHH(theta, model)
  meanTotal <- mu$mean

  ############################################################
  
  ## helper function for derivatives
  derivHHH.factor <- if(dimPsi > 0L){ # NegBin
      psiPlusMu <- psi + meanTotal    # also used below for calculation of grPsi
      psiYpsiMu <- (psi+Y) / psiPlusMu
      Y/meanTotal - psiYpsiMu
  } else { # Poisson
      Y/meanTotal - 1
  }
  derivHHH <- function (dmu) derivHHH.factor * dmu

  ## go through groups of parameters and compute the gradient of each component
  computeGrad <- function(mean.comp){
  
    grad.fe <- numeric(0L)
    grad.re <- numeric(0L)
       
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
        grad.re <- c(grad.re, dRTheta)
        grad.fe <- c(grad.fe, sum(dTheta))
      } else{
        grad.fe <- c(grad.fe, summ(dTheta,na.rm=TRUE))
      }
    }
    
    list(fe=grad.fe, re=unname(grad.re))
  }
  
  gradients <- computeGrad(mu[c("epi.own","epi.neighbours","endemic")])

  ## gradient for parameter vector of the neighbourhood weights
  grd <- if (dimd > 0L) {
      dneOffset <- model$offset[[2L]](pars$d, type="gradient")
      ##<- this is always a list (of length dimd) of matrices/arrays
      onescore.d <- function (dneoff) {
          dmudd <- mu$ne.exppred * dneoff[subset,,drop=FALSE]
          grd.terms <- derivHHH(dmudd)
          sum(grd.terms, na.rm=TRUE)
      }
      unlist(clapply(dneOffset, onescore.d), recursive=FALSE, use.names=FALSE)
  } else numeric(0L)
  
  ## gradient for overdispersion parameter psi
  grPsi <- if(dimPsi > 0L){
      dPsi <- psi * (digamma(Y+psi) - digamma(psi) + log(psi) + 1
                     - log(psiPlusMu) - psiYpsiMu)
      # multiple psi_i's or single psi?
      if(dimPsi>1L) colSums(dPsi, na.rm=TRUE) else sum(dPsi, na.rm=TRUE)
  } else numeric(0L)
  
  if(any(is.na(grPsi))){
      stop("derivatives for overdispersion parameter psi not computable")
  }
  
  ## add penalty to random effects gradient
  s.pen <- if(dimRE > 0) c(Sigma.inv %*% randomEffects) else numeric(0L)
  if(length(gradients$re) != length(s.pen)) # FIXME: could this ever happen?
      cat("length of s(b) != Sigma.inv %*% b\n") # should probably better stop()
  grRandom <- c(gradients$re - s.pen)

  ## Done
  res <- c(gradients$fe, grd, grPsi, grRandom)
  res
}


penFisher <- function(theta, sd.corr, model, attributes=FALSE)
{
  if(any(is.na(theta))) stop("NAs in regression parameters.", ADVICEONERROR)

  ## unpack model
  subset <- model$subset
  Y <- model$response[subset,,drop=FALSE]
  isNA <- model$isNA[subset,,drop=FALSE]
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE

  ## unpack parameters
  pars <- splitParams(theta, model)
  psi <- exp(pars$overdisp)             # = 1/psi, pars$overdisp = -log(psi)
  if(dimPsi > 1L) {
      psi <- matrix(psi, ncol=model$nUnits, nrow=model$nTime,
                    byrow=TRUE)[subset,,drop=FALSE]
  }
  if (dimRE > 0) {
      randomEffects <- pars$random
      sd   <- head(sd.corr, model$nVar)
      corr <- tail(sd.corr, model$nCorr)
      dimBlock <- model$rankRE[model$rankRE>0]
      Sigma.inv <- getSigmaInv(sd, corr, model$nVar, dimBlock)
  }

  ## further model unpacking
  term <- model$terms
  nGroups <- model$nGroups
  dimd  <- model$nd
  dimFE <- model$nFE
  idxFE <- model$indexFE
  idxRE <- model$indexRE

  ## evaluate mean
  mu <- meanHHH(theta, model)
  meanTotal <- mu$mean
  
  ############################################################
  
  ## helper functions for derivatives:
  if (dimPsi > 0L) { # negbin
    psiPlusY <- psi + Y
    psiPlusMu <- psi + meanTotal
    psiPlusMu2 <- psiPlusMu^2
    psiYpsiMu <- psiPlusY / psiPlusMu
    psiYpsiMu2 <- psiPlusY / psiPlusMu2
    deriv2HHH.fac1 <- psiYpsiMu2 - Y / (meanTotal^2)
    deriv2HHH.fac2 <- Y / meanTotal - psiYpsiMu
    ## psi-related derivatives
    dThetadPsi.fac <- psi * (psiYpsiMu2 - 1/psiPlusMu)
    dThetadPsi <- function(dTheta){
        dThetadPsi.fac * dTheta
    }
    dPsidPsi <- function(){
        dPsi <- psi * (digamma(psiPlusY)-digamma(psi) +log(psi)+1 -
                       log(psiPlusMu) - psiYpsiMu)
        psi^2 * (trigamma(psiPlusY) - trigamma(psi) + 1/psi - 1/psiPlusMu -
                 (meanTotal-Y)/psiPlusMu2) + dPsi
    }
  } else { # poisson
    deriv2HHH.fac1 <- -Y / (meanTotal^2)
    deriv2HHH.fac2 <- Y / meanTotal - 1
  }
  deriv2HHH <- function(dTheta_l, dTheta_k, dTheta_lk){
      dTheta_l * dTheta_k * deriv2HHH.fac1 + dTheta_lk * deriv2HHH.fac2
  }

  ## go through groups of parameters and compute the hessian of each component
  computeFisher <- function(mean.comp){
    # initialize hessian
    hessian.FE.FE <- matrix(0,dimFE,dimFE)
    hessian.FE.RE <- matrix(0,dimFE,dimRE)
    hessian.RE.RE <- matrix(0,dimRE,dimRE)
    
    hessian.FE.Psi <- matrix(0,dimFE,dimPsi)
    hessian.Psi.RE <- matrix(0,dimPsi,dimPsi+dimRE) # CAVE: contains PsiPsi and PsiRE

    hessian.FE.d <- matrix(0,dimFE,dimd)
    hessian.d.d <- matrix(0,dimd,dimd)
    hessian.d.Psi <- matrix(0,dimd,dimPsi)
    hessian.d.RE <- matrix(0,dimd,dimRE)
    
    ## derivatives wrt neighbourhood weight parameters d
    if (dimd > 0L) {
        phi.doff <- function (dneoff) {
            mu$ne.exppred * dneoff[subset,,drop=FALSE]
        }
        ## for type %in% c("gradient", "hessian"), model$offset[[2L]] is always
        ## a list of matrices/arrays. It has length(pars$d) elements for the
        ## gradient and length(pars$d)*(length(pars$d)+1)/2 for the hessian.
        dneOffset <- model$offset[[2L]](pars$d, type="gradient")
        dmudd <- lapply(dneOffset, phi.doff)
        d2neOffset <- model$offset[[2L]](pars$d, type="hessian")
        d2mudddd <- lapply(d2neOffset, phi.doff)
        ## d l(theta,x) /dd dd (fill only upper triangle)
        ij <- 0L
        for (i in seq_len(dimd)) {
            for (j in i:dimd) {
                ij <- ij + 1L  #= dimd*(i-1) + j - (i-1)*i/2  # for j >= i
                ## d2mudddd contains upper triangle by row (=lowertri by column)
                d2ij <- deriv2HHH(dmudd[[i]], dmudd[[j]], d2mudddd[[ij]])
                hessian.d.d[i,j] <- sum(d2ij, na.rm=TRUE)
            }
        }
    }

    if (dimPsi > 0L) {
        ## d l(theta,x) /dpsi dpsi
        hessian.Psi.RE[,seq_len(dimPsi)] <- if (dimPsi == 1L) {
            sum(dPsidPsi(), na.rm=TRUE)
        } else diag(colSums(dPsidPsi(), na.rm=TRUE))
        ## d l(theta) / dd dpsi
        for (i in seq_len(dimd)) {      # will not be run if dimd==0
            dPsi.i <- colSums(dThetadPsi(dmudd[[i]]),na.rm=TRUE)
            if(dimPsi==1L) dPsi.i <- sum(dPsi.i)
            hessian.d.Psi[i,] <- dPsi.i
        }
    }
    
    ##
    i.fixed <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]  
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums(didj %mj% Z.j)
        ##<- didj must not contain NA's (all NA's set to 0)
        dIJ <- sum(didj,na.rm=TRUE)     # fixed on 24/09/2012
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
        "%mj%" <- get(term["mult",j][[1]]) 
        dIJ <- colSums(didj %mj% Z.j)   # didj must not contain NA's (all NA's set to 0)
        hessian.FE.RE[idxFE==i,idxRE==j] <<- diag(dIJ)[ which.i, ]
        dIJ <- dIJ[ which.i ]           # added which.i subsetting in r432
                                        # FIXME @Michaela: please confirm this correction
      } else if(unitSpecific.j){
        which.ij <- (which.i & which.j) # FIXME: this is actually unused...?
        dIJ <- diag(colSums(didj))[ which.i, which.j ] 
      } else {
        dIJ <- colSums(didj)[ which.i ]
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ  
    }
    ##
    i.random <- function(){
      if(random.j){
        Z.j <- term["Z.intercept",j][[1]]  
        "%mj%" <- get(term["mult",j][[1]])
        hessian.FE.RE[idxFE==i,idxRE==j] <<- colSums(didj %mj% Z.j)
        hessian.FE.RE[idxFE==j,idxRE==i] <<- colSums(didj %m% Z.i)

        if(length(Z.j)==1 & length(Z.i)==1){
          Z <- Z.i*Z.j
          hessian.RE.RE[which(idxRE==i),idxRE==j] <<- diag(colSums( didj %m% Z))
        } else if(length(Z.j)==1 & length(Z.i)>1){         #*
          Z.j <- diag(nrow=model$nUnits)
          for(k in 1:ncol(Z.j)){
            Z <- Z.i*Z.j[,k]
            hessian.RE.RE[idxRE==i,which(idxRE==j)[k]] <<- colSums( didj %m% Z)
          }  
        } else if(length(Z.j)>1 & length(Z.i)==1){         #* 
          Z.i <- diag(nrow=model$nUnits)
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
        dIJ <- sum(didj)
      } else if(unitSpecific.j){
        dIJ <- colSums(didj %m% Z.i)
        hessian.FE.RE[idxFE==j,idxRE==i] <<- diag(dIJ)[ which.j, ]
        dIJ <- dIJ[ which.j ]
      } else {
        hessian.FE.RE[idxFE==j,idxRE==i] <<- colSums(didj %m% Z.i)
        dIJ <- sum(didj)
      }
      hessian.FE.FE[idxFE==i,idxFE==j] <<- dIJ    
    }
  ##----------------------------------------------           

    for(i in 1:nGroups){ #go through rows of hessian
      # parameter group belongs to which components
      comp.i <- term["offsetComp",i][[1]]      
      # get covariate value
      Xit <- term["terms",i][[1]] # eiter 1 or a matrix with values
      if(is.matrix(Xit)){
        Xit <- Xit[subset,,drop=FALSE]
      }
      m.Xit <- mean.comp[[comp.i]]*Xit
      
      random.i <- term["random",i][[1]]
      unitSpecific.i <- term["unitSpecific",i][[1]]

      ## fill psi-related entries and select fillHess function
      if(random.i){
        Z.i <- term["Z.intercept",i][[1]]   # Z.i and %m% (of i) determined here
        "%m%" <- get(term["mult",i][[1]])   # will also be used in j's for loop
        fillHess <- i.random
        if(dimPsi==1L){
          dPsi<- dThetadPsi(m.Xit)
          dPsi[isNA] <- 0
          hessian.FE.Psi[idxFE==i,]<- sum(dPsi)
          hessian.Psi.RE[,c(FALSE,idxRE==i)] <- colSums( dPsi%m%Z.i)
        } else if(dimPsi>1L){
          dPsi<- dThetadPsi(m.Xit)
          dPsi[isNA] <- 0
          hessian.FE.Psi[idxFE==i,]<- colSums( dPsi)
          hessian.Psi.RE[,c(rep.int(FALSE, dimPsi),idxRE==i)] <- diag(colSums( dPsi%m%Z.i))
        }
        
      } else if(unitSpecific.i){
        which.i <- term["which",i][[1]]
        fillHess <- i.unit
        if(dimPsi>0L){
          dPsi <- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
          hessian.FE.Psi[idxFE==i,] <- if(dimPsi==1L) {
              dPsi[which.i]
          } else {
              diag(dPsi)[which.i,]
          }
        }
        
      } else {
        fillHess <- i.fixed
        if(dimPsi>0L){
          dPsi<- colSums(dThetadPsi(m.Xit),na.rm=TRUE)
          if(dimPsi==1L){ dPsi <- sum(dPsi) }
          hessian.FE.Psi[idxFE==i,] <-dPsi
        }
      }

      ## fill pars$d-related entries
      for (j in seq_len(dimd)) {      # will not be run if dimd==0
          didd <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = dmudd[[j]],
                            dTheta_lk = if (comp.i == 2) dmudd[[j]] * Xit else 0)
          didd[isNA] <- 0
          hessian.FE.d[idxFE==i,j] <- if (unitSpecific.i) {
              colSums(didd,na.rm=TRUE)[which.i]
          } else sum(didd)
          if (random.i) hessian.d.RE[j,idxRE==i] <- colSums(didd %m% Z.i)
      }
      
      ## fill other (non-psi, non-d) entries (only upper triangle, j >= i!)
      for(j in i:nGroups){
        comp.j <- term["offsetComp",j][[1]]
        
        Xjt <- term["terms",j][[1]] # eiter 1 or a matrix with values
        if(is.matrix(Xjt)){
          Xjt <- Xjt[subset,,drop=FALSE]
        }
        # if param i and j do not belong to the same component, d(i)d(j)=0
        m.Xit.Xjt <- if (comp.i != comp.j) 0 else m.Xit * Xjt
        
        didj <- deriv2HHH(dTheta_l = m.Xit, dTheta_k = mean.comp[[comp.j]]*Xjt,
                          dTheta_lk = m.Xit.Xjt)
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
    hessian <- rbind(cbind(hessian.FE.FE,hessian.FE.d,hessian.FE.Psi,hessian.FE.RE),
                     cbind(matrix(0,dimd,dimFE),hessian.d.d,hessian.d.Psi,hessian.d.RE),
                     cbind(matrix(0,dimPsi,dimFE+dimd),hessian.Psi.RE),
                     cbind(matrix(0,dimRE,dimFE+dimd+dimPsi),hessian.RE.RE))

    hessian[lower.tri(hessian)] <- 0    # FIXME: should already be the case!
    diagHessian <- diag(hessian)
    fisher <- -(hessian + t(hessian))
    diag(fisher) <- -diagHessian   
  
    return(fisher)
  }

  fisher <- computeFisher(mu[c("epi.own","epi.neighbours","endemic")])

  ## add penalty for random effects
  pen <- matrix(0, length(theta), length(theta))
  Fpen <- if(dimRE > 0){
    thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
    pen[thetaIdxRE,thetaIdxRE] <- Sigma.inv
    fisher + pen
  } else fisher

  ## Done
  if(attributes){
    attr(Fpen, "fisher") <- fisher
    attr(Fpen, "pen") <- pen
  }
  Fpen
}


#################################################

sqrtOf1pr2 <- function(r){ sqrt(1+r^2) }

getSigmai <- function(sd,    # vector of length dim with stdev's
                      correlation, # vector of length dim with correlation parameters, =0 if un-correlated
                      dim
                      ){
  if(dim==0) return(NULL)
                      
  Sigma.i <- if (length(correlation) == 0L) diag(exp(2*sd), dim) else {
      D <- diag(exp(sd), dim)
      L <- diag(nrow=dim)
      L[2,1:2] <- c(correlation[1],1)/sqrtOf1pr2(correlation[1])
      if (dim==3) {
          L[3,] <- c(correlation[2:3],1)/sqrtOf1pr2(correlation[2])
          L[3,2:3] <- L[3,2:3]/sqrtOf1pr2(correlation[3])
      }
      D %*% tcrossprod(L) %*% D  # ~75% quicker than D %*% L %*% t(L) %*% D
  }
  return(Sigma.i)
}

getSigmaiInv <- function(sd,    # vector of length dim with stdev's
                         correlation, # vector of length dim with correlation parameters, =0 if un-correlated
                         dim
                         ){
    
  if(dim==0) return(NULL) 
  
  Sigma.i.inv <- if (length(correlation) == 0L) diag(exp(-2*sd), dim) else {
      r <- correlation
      Dinv <- diag(exp(-sd), dim)
      L <- diag(nrow=dim)
      L[2,1:2] <- c(-r[1],sqrtOf1pr2(r[1]))
      if(dim==3){
          L[3,1] <- r[1]*r[3]-r[2]*sqrtOf1pr2(r[3])
          L[3,2] <- -L[2,2]*r[3]
          L[3,3] <- sqrtOf1pr2(r[2])*sqrtOf1pr2(r[3])
      }
      Dinv %*% crossprod(L) %*% Dinv  # ~75% quicker than Dinv %*% t(L) %*% L %*% Dinv
  }
  
  return(Sigma.i.inv)
}

#* allow blockdiagonal matrix blockdiag(A,B), with A=kronecker product, B=diagonal matrix?
getSigmaInv <- function(sd, correlation, dimSigma, dimBlocks, SigmaInvi=NULL){
  if(is.null(SigmaInvi)){
    SigmaInvi <- getSigmaiInv(sd,correlation,dimSigma)
  }
  SigmaInv <- if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    kronecker(SigmaInvi,diag(nrow=dimBlocks[1]))
    # the result is a symmetric matrix if SigmaInvi is symmetric
  } else { # kronecker product not possible -> correlation=0 
    diag(rep(diag(SigmaInvi),dimBlocks))
  }
  return(SigmaInv)
}

getSigma <- function(sd, correlation, dimSigma, dimBlocks, Sigmai=NULL){
  if(is.null(Sigmai)){
    Sigmai <- getSigmai(sd,correlation,dimSigma)
  }
  Sigma <- if(length(unique(dimBlocks))==1){  # kronecker product formulation possible
    kronecker(Sigmai,diag(nrow=dimBlocks[1]))
    # the result is a symmetric matrix if Sigmai is symmetric
  } else { # kronecker product not possible -> correlation=0 
    diag(rep(diag(Sigmai),dimBlocks))
  }
  return(Sigma)
}


## Approximate marginal likelihood for variance components
## Parameter and model unpacking at the beginning (up to the ###...-line) is
## identical in marScore() and marFisher()
marLogLik <- function(sd.corr, theta, model, fisher.unpen=NULL){

  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(-Inf)
  }
  
  if(any(is.na(sd.corr))){
    # in order to avoid nlminb from running into an infinite loop (cf. bug
    # report #15052), we have to emergency stop() in this case.
    # As of R 2.16.0, nlminb() will throw an error if it receives NA from
    # any of the supplied functions.
    stop("NAs in variance parameters.", ADVICEONERROR)
  }

  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- splitParams(theta,model)  
  randomEffects <- pars$random

  dimFE <- model$nFE
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  
  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  ############################################################
  
  # penalized part of likelihood
  # compute -0.5*log(|Sigma|) - 0.5*RE' %*% Sigma.inv %*% RE
  # where -0.5*log(|Sigma|) = -dim(RE_i)*[Sum(sd_i) -0.5*log(1+corr_i^2)]
  ##lpen <- -0.5*(t(randomEffects)%*%Sigma.inv%*%randomEffects)
  ## the following implementation takes ~85% less computing time !
  lpen <- -0.5 * crossprod(randomEffects, Sigma.inv) %*% randomEffects
  loglik.pen <- sum(-dimBlocks*sd) + c(lpen)
  if(dimCorr >0){
    loglik.pen <- loglik.pen + 0.5*dimBlocks[1]*sum(log(1+corr^2))
  }
  
  ## approximate marginal likelihood
  logdetfisher <- determinant(fisher,logarithm=TRUE)$modulus
  lmarg <- loglik.pen -0.5*c(logdetfisher)
      
  return(lmarg)
}


marScore <- function(sd.corr, theta,  model, fisher.unpen=NULL){

  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(numeric(0L))
  }

  if(any(is.na(sd.corr))) stop("NAs in variance parameters.", ADVICEONERROR)
  
  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- splitParams(theta,model)  
  randomEffects <- pars$random

  dimFE <- model$nFE
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  
  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  # inverse of penalized fisher info
  F.inv <- try(solve(fisher),silent=TRUE)
  if(inherits(F.inv,"try-error")){
    cat("\n WARNING (in marScore): penalized Fisher is singular!\n")
    #return(rep.int(0,dimSigma))
    ## continuing with the generalized inverse often works, otherwise we would
    ## have to stop() here, because nlminb() cannot deal with NA's
    F.inv <- MASS::ginv(fisher)
  }
  F.inv.RE <- F.inv[thetaIdxRE,thetaIdxRE]

  ############################################################
  
  ## compute marginal score and fisher for each variance component
  # initialize score and fisher info
  marg.score <- rep.int(NA_real_,dimSigma)
  
  ## specify functions for derivatives
  deriv1 <- switch(dimVar, dSigma1, dSigma2, dSigma3)
  
  d1Sigma <- deriv1(sd, corr)
  Sigmai.inv <- getSigmaiInv(sd, corr, dimVar)
  
 
  # derivation of log determinant
  # -.5*tr(Sigma^-1 %*% dSigma/ds) = -R (for sd.i) 
  #                                =  R*corr.i/(corr.i^2+1) (for corr.i)
  d1logDet <- c(-dimBlocks,dimBlocks[1]*corr/(corr^2+1))
  
  # go through all variance parameters
  for(i in 1:dimSigma){
    dSi <- -Sigmai.inv %*% d1Sigma[,,i] %*% Sigmai.inv # CAVE: sign
    dS.i <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSi)
    #dlpen.i <- -0.5* t(randomEffects) %*% dS.i %*% randomEffects
    # ~85% faster implementation using crossprod() avoiding "slow" t():
    dlpen.i <- -0.5 * crossprod(randomEffects, dS.i) %*% randomEffects
    #tr.d1logDetF <- sum(diag(F.inv.RE %*% dS.i))
    tr.d1logDetF <- sum(F.inv.RE * dS.i)   # since dS.i is symmetric
    #<- needs 1/100 (!) of the computation time of sum(diag(F.inv.RE %*% dS.i))
    marg.score[i] <- d1logDet[i] + c(dlpen.i) - 0.5 * tr.d1logDetF
  }
  
  return(marg.score)
}


marFisher <- function(sd.corr, theta,  model, fisher.unpen=NULL){
  
  dimVar <- model$nVar
  dimCorr <- model$nCorr
  dimSigma <- model$nSigma
  
  if(dimSigma == 0){
    return(matrix(numeric(0L),0L,0L))
  }

  if(any(is.na(sd.corr))) stop("NAs in variance parameters.", ADVICEONERROR)
  
  sd   <- head(sd.corr,dimVar)
  corr <- tail(sd.corr,dimCorr)

  pars <- splitParams(theta,model)  
  randomEffects <- pars$random

  dimFE <- model$nFE
  dimPsi <- model$nOverdisp
  dimRE <- model$nRE
  
  dimBlocks <- model$rankRE[model$rankRE>0]
  Sigma.inv <- getSigmaInv(sd, corr, dimVar, dimBlocks)

  # if not given, calculate unpenalized part of fisher info 
  if(is.null(fisher.unpen)){
    fisher.unpen <- attr(penFisher(theta, sd.corr, model,attributes=TRUE), "fisher")
  }
  
  # add penalty to fisher
  fisher <- fisher.unpen
  thetaIdxRE <- seq.int(to=length(theta), length.out=dimRE)
  fisher[thetaIdxRE,thetaIdxRE] <- fisher[thetaIdxRE,thetaIdxRE] + Sigma.inv

  # inverse of penalized fisher info
  F.inv <- try(solve(fisher),silent=TRUE)
  if(inherits(F.inv,"try-error")){
    cat("\n WARNING (in marFisher): penalized Fisher is singular!\n")
    #return(matrix(Inf,dimSigma,dimSigma))
    ## continuing with the generalized inverse often works, otherwise we would
    ## have to stop() here, because nlminb() cannot deal with NA's
    F.inv <- MASS::ginv(fisher)
  }
  F.inv.RE <- F.inv[thetaIdxRE,thetaIdxRE]
  ## declare F.inv.RE as a symmetric matrix?
  ##F.inv.RE <- new("dsyMatrix", Dim = dim(F.inv.RE), x = c(F.inv.RE))
  ## -> no, F.inv.RE %*% dS.i becomes actually slower (dS.i is a "sparseMatrix")
  
  ############################################################
  
  marg.hesse <- matrix(NA_real_,dimSigma,dimSigma)
  
  ## specify functions for derivatives
  deriv1 <- switch(dimVar,dSigma1, dSigma2, dSigma3)
  deriv2 <- switch(dimVar,d2Sigma1, d2Sigma2, d2Sigma3)
  
  d1Sigma <- deriv1(sd, corr)
  d2Sigma <- deriv2(sd, corr, d1Sigma)
  Sigmai.inv <- getSigmaiInv(sd, corr, dimVar)
  
 
  # 2nd derivatives of log determinant
  d2logDet <- diag(c(rep.int(0,dimVar),-dimBlocks[1]*(corr^2-1)/(corr^2+1)^2),dimSigma)

  # function to convert dS.i and dS.j matrices to sparse matrix objects
  dS2sparse <- if (dimCorr > 0) function (x) {
      Matrix::forceSymmetric(as(x, "sparseMatrix")) # dS.i & dS.j are symmetric
  } else function (x) { #as(x, "diagonalMatrix")
      new("ddiMatrix", Dim = dim(x), diag = "N", x = diag(x))
  }
  
  # go through all variance parameters
  for(i in 1:dimSigma){
    # compute first derivative of the penalized Fisher info (-> of Sigma^-1) 
    # with respect to the i-th element of Sigma (= kronecker prod. of Sigmai and identity matrix)
    # Harville Ch15, Eq. 8.15: (d/d i)S^-1 = - S^-1 * (d/d i) S * S^-1
    SigmaiInv.d1i <- Sigmai.inv %*% d1Sigma[,,i]
    dSi <- -SigmaiInv.d1i %*% Sigmai.inv
    dS.i <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSi)
    dS.i <- dS2sparse(dS.i)
    
    # compute second derivatives
    for(j in i:dimSigma){
      # compute (d/d j) S^-1
      SigmaiInv.d1j <- Sigmai.inv %*% d1Sigma[,,j]
      dSj <- -SigmaiInv.d1j %*% Sigmai.inv
      dS.j <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=dSj)
      dS.j <- dS2sparse(dS.j)
      # compute (d/di dj) S^-1
      #dS.ij <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,
      #                  Sigmai=d2Sigma[[i]][,,j])
      
      # compute second derivatives of Sigma^-1 (Harville Ch15, Eq 9.2)
      d2S <- (- Sigmai.inv %*% d2Sigma[[i]][,,j] + 
               SigmaiInv.d1i %*% SigmaiInv.d1j +
               SigmaiInv.d1j %*% SigmaiInv.d1i) %*% Sigmai.inv 
      
      dSij <- getSigma(dimSigma=dimVar,dimBlocks=dimBlocks,Sigmai=d2S)

      #d2lpen.i <- -0.5* t(randomEffects) %*% dSij %*% randomEffects
      # ~85% faster implementation using crossprod() avoiding "slow" t():
      d2lpen.i <- -0.5 * crossprod(randomEffects, dSij) %*% randomEffects

      mpart1 <- dS.j %*% F.inv.RE  # 3 times as fast as the other way round
      mpart2 <- dS.i %*% F.inv.RE
      mpart <- mpart1 %*% mpart2
      ## speed-ups: - tr(F.inv.RE %*% dSij) simply equals sum(F.inv.RE * dSij)
      ##            - accelerate matrix product by sparse matrices dS.i and dS.j
      ##            - use cyclic permutation of trace:
      ##              tr(F.inv.RE %*% dS.j %*% F.inv.RE %*% dS.i) =
      ##              tr(dS.j %*% F.inv.RE %*% dS.i %*% F.inv.RE)
      tr.d2logDetF <- -sum(Matrix::diag(mpart)) + sum(F.inv.RE * dSij)
      
      marg.hesse[i,j] <- marg.hesse[j,i] <-
          d2logDet[i,j] + d2lpen.i - 0.5 * tr.d2logDetF
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
    result$dcorr1[2,1,3] <- result$dcorr1[1,2,3] <-
        -(3*corr[1]*exp(sum(sd[1:2])))/(sqrtOf1pr2(corr[1])^5)
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


updateRegression <- function(theta, sd.corr, model, control, method = "nlminb")
{  
  lower <- control[["lower"]]; control$lower <- NULL
  upper <- control[["upper"]]; control$upper <- NULL
  if(any(theta %in% lower[is.finite(lower)])){
    cat("WARNING: at least one theta reached lower bound\n")
  }
  
  # estimate regression coefficients (theta, fixed + random)
  if(method == "nlminb"){
    ll <- function(theta,...) {- penLogLik(theta,...)}
    gr <- function(theta,...) {- penScore(theta,...)}
    scale <- control[["scale"]]; control$scale <- NULL
    res <- nlminb(start=theta, objective=ll, gradient=gr, hessian=penFisher, 
                  sd.corr=sd.corr, model=model,
                  scale=scale, control=control, lower=lower, upper=upper)
    theta.new <- res$par
    ll <- - res$objective
  } else if(method=="nr"){
    # objective function
    regressionParams <- function(theta,...){
      ll <- penLogLik(theta,...)
      attr(ll,"score") <- penScore(theta,...)
      attr(ll,"fisher") <- penFisher(theta,...)
      return(ll)
    }
    res <- newtonRaphson(x=theta, fn=regressionParams,
                         sd.corr=sd.corr, model=model,
                         control=control, verbose=control$verbose)
    theta.new <- res$coefficients
    ll <- res$loglik
  } else if(method == "nlm"){
    # objective function
    regressionParamsMin <- function(theta,...){
      ll <- - penLogLik(theta,...)
      attr(ll,"gradient") <- - penScore(theta,...)
      attr(ll,"hessian") <- penFisher(theta,...)
      return(ll)  
    }
    res <- do.call("nlm", args=c(alist(p=theta, f=regressionParamsMin,
                          sd.corr=sd.corr, model=model), control))
    theta.new <- res$estimate
    ll <- - res$minimum
    # nlm returns convergence status in $code, 1-2 indicate convergence,
    # 3-5 indicate non-convergence
    res$convergence <- as.numeric(res$code >2)
  } else { # use optim
    res <- optim(theta, penLogLik, penScore, sd.corr=sd.corr, model=model,
                 method=method, lower=lower, upper=upper, control=control)
    theta.new <- res$par
    ll <- res$value
  }

  rel.tol <- max(abs(theta.new - theta)) / max(abs(theta))
  # The above is a weaker criterion than the maximum relative parameter change:
  #rel.tol <- max(abs(theta.new/theta - 1))
  
  return(list(par=theta.new, ll=ll, rel.tol=rel.tol,
              convergence=res$convergence, message=res$message))
}


updateVariance <- function(sd.corr, theta, model, fisher.unpen,
                           control, method = "nlminb")
{
  # only fixed effects => no variance 
  if(length(sd.corr) == 0L){
    return(list(par=sd.corr,ll=NA_real_,rel.tol=0,convergence=0))
  }

  lower <- control[["lower"]]; control$lower <- NULL
  upper <- control[["upper"]]; control$upper <- NULL
  
  # estimate variance parameters (sd.corr)
  if(method == "nlminb"){
    ll <- function(sd.corr,...) { - marLogLik(sd.corr,...)}
    gr <- function(sd.corr,...) { - marScore(sd.corr,...)}
    scale <- control[["scale"]]; control$scale <- NULL
    res <- nlminb(start=sd.corr, objective=ll, gradient=gr, hessian=marFisher,
                  theta=theta, model=model, fisher.unpen=fisher.unpen,
                  scale=scale, control=control, lower=lower, upper=upper)
    sd.corr.new <- res$par
    ll <- -res$objective
    if(any(is.na(sd.corr.new))) {
        ## Before R 2.16.0, in some cases nlminb continues being stucked at the
        ## NA parameters until the iteration limit is reached (then returning).
        ## In others, nlminb is in an infinite loop in the FORTRAN code even
        ## ignoring the BREAK signal.
        if (res$convergence == 0) res$convergence <- 99
    }
  } else if(method == "nr"){
    varianceParams <-  function(sd.corr,...){
      ll <- marLogLik(sd.corr,...)
      attr(ll,"score") <- marScore(sd.corr,...)
      attr(ll,"fisher") <- marFisher(sd.corr,...)
      return(ll)
    }
    res <- newtonRaphson(x=sd.corr, fn=varianceParams,
                         theta=theta, model=model, fisher.unpen=fisher.unpen,
                         control=control, verbose=control$verbose)
    sd.corr.new <- res$coefficients
    ll <- res$loglik
  } else if(method == "nlm"){
    varianceParamsMin <-  function(sd.corr,...){
      ll <- - marLogLik(sd.corr,...)
      attr(ll,"gradient") <- - marScore(sd.corr,...)
      attr(ll,"hessian") <- marFisher(sd.corr,...)
      return(ll)
    }
    res <- do.call("nlm", args=c(alist(p=sd.corr, f=varianceParamsMin,
                          theta=theta, model=model, fisher.unpen=fisher.unpen),
                          control))
    sd.corr.new <- res$estimate
    ll <- -res$minimum
    # nlm returns convergence status in $code, 1-2 indicate convergence,
    # 3-5 indicate non-convergence
    res$convergence <- as.numeric(res$code >2)
  } else { # use optim
    res <- optim(sd.corr, marLogLik, marScore,
                 theta=theta, model=model, fisher.unpen=fisher.unpen,
                 method=method, lower=lower, upper=upper, control=control)
    sd.corr.new <- res$par  
    ll <- res$value
  }

  rel.tol <- max(abs(sd.corr.new - sd.corr)) / max(abs(sd.corr)) 
  # The above is a weaker criterion than the maximum relative parameter change:
  #rel.tol <- max(abs(sd.corr.new/sd.corr - 1))
  
  return(list(par=sd.corr.new, ll=ll, rel.tol=rel.tol,
              convergence=res$convergence, message=res$message))
}


## default control arguments for updates
defaultOptimControl <- function (method = "nlminb", lower = -Inf, upper = Inf,
                                 verbose = 0)
{
    lowVerbose <- verbose %in% 0:2
    luOptimMethod <- method %in% c("Brent", "L-BFGS-B")
    defaults.nr <- list(scoreTol=1e-5, paramTol=1e-7, F.inc=0.01, stepFrac=0.5,
                        niter=20, verbose=verbose)
    defaults.nlminb <- list(iter.max=20, scale=1, lower=lower, upper=upper,
                            trace=if(lowVerbose) c(0,0,5)[verbose+1] else 1)
    defaults.nlm <- list(iterlim=100, check.analyticals=FALSE,
                         print.level=if(lowVerbose) c(0,0,1)[verbose+1] else 2)
    defaults.optim <- list(fnscale=-1, trace=max(0,verbose-1),
                           lower=if (luOptimMethod) lower else -Inf,
                           upper=if (luOptimMethod) upper else Inf)
    defaults <- switch(method, "nr" = defaults.nr, "nlm" = defaults.nlm,
                       "nlminb" = defaults.nlminb, defaults.optim)
    return(defaults)
}


## fitHHH is the main workhorse where the iterative optimization is performed
fitHHH <- function(theta, sd.corr, model,
                   cntrl.stop=list(tol=1e-5, niter=100),
                   cntrl.regression=list(method="nlminb"),
                   cntrl.variance=list(method="nlminb"),
                   verbose=0, shrinkage=FALSE)
{  
  dimFE.d.O <- model$nFE + model$nd + model$nOverdisp
  dimRE <- model$nRE

  ## artificial lower bound on regression parameters
  reg.lower <- rep.int(-Inf, length(theta))
  indexAR <- c(grep("ar.ri",model$namesFE), grep("ar.1",model$namesFE),
               grep("ne.ri",model$namesFE), grep("ne.1",model$namesFE))
  reg.lower[indexAR] <- -20

  ## control arguments
  reg.method <- cntrl.regression$method; cntrl.regression$method <- NULL
  var.method <- cntrl.variance$method; cntrl.variance$method <- NULL
  if (length(sd.corr) == 1L && var.method == "Nelder-Mead") {
      var.method <- "Brent"
      cat("NOTE: switched variance optimizer from", dQuote("Nelder-Mead"),
           "to", dQuote("Brent"), "(dim(Sigma)=1)\n")
  }
  defaults.reg <- defaultOptimControl(method=reg.method, lower=reg.lower,
                                      upper=Inf, verbose=verbose)
  defaults.var <- defaultOptimControl(method=var.method, lower=-5, upper=5,
                                      verbose=verbose)
  cntrl.update.reg <- modifyList(defaults.reg, cntrl.regression)
  cntrl.update.var <- modifyList(defaults.var, cntrl.variance)
  if (!is.null(cntrl.update.reg$fnscale) && cntrl.update.reg$fnscale > 0) {
      cntrl.update.reg$fnscale <- -cntrl.update.reg$fnscale
  }
  if (!is.null(cntrl.update.var$fnscale) && cntrl.update.var$fnscale > 0) {
      cntrl.update.var$fnscale <- -cntrl.update.var$fnscale
  }

  ## Let's go
  if (verbose>0) {
      cat(as.character(Sys.time()), ":",
          "Iterative optimization of regression & variance parameters\n")
  }
  convergence <- 99
  i <- 0
  while(convergence != 0 && (i < cntrl.stop$niter)){
    i <- i+1
    if (verbose>0) cat("\n")
    
    ## update regression coefficients
    parReg <- updateRegression(theta=theta, sd.corr=sd.corr, model=model,
                               control=cntrl.update.reg, method=reg.method)
    theta <- parReg$par
    fisher.unpen <- attr(penFisher(theta, sd.corr, model, attributes=TRUE),
                         "fisher")

    if(verbose>0)
      cat("Update of regression parameters:  max|x_0 - x_1| / max|x_0| =",
          parReg$rel.tol, "\n")
    
    if(parReg$convergence!=0) {
      if (!is.null(parReg$message)) print(parReg$message)
      cat("Update of regression coefficients in iteration ", i, " unreliable\n")
    }

    if(parReg$convergence >20 && shrinkage){
      cat("\n\n***************************************\nshrinkage",
          0.1*theta[abs(theta)>10],"\n")
      theta[abs(theta)>10] <- 0.1*theta[abs(theta)>10]
      diag(fisher.unpen) <- diag(fisher.unpen)+1e-2
    }
    
    ## update variance parameters
    parVar <- updateVariance(sd.corr=sd.corr, theta=theta, model=model,
                             fisher.unpen=fisher.unpen,
                             control=cntrl.update.var, method=var.method)

    if(verbose>0)
      cat("Update of variance parameters:  max|x_0 - x_1| / max|x_0| =",
          parVar$rel.tol, "\n")
      
    if(parVar$convergence!=0) {
      if (!is.null(parVar$message)) print(parVar$message)
      cat("Update of variance parameters in iteration ", i, " unreliable\n")
    }

    if(any(is.na(parVar$par))){         # this might occur when using nlminb
      var.method <- if (length(sd.corr) == 1L) "Brent" else "Nelder-Mead"
      cat("WARNING: at least one updated variance parameter is not a number\n",
          "\t-> NO UPDATE of variance\n",
          "\t-> SWITCHING to robust", dQuote(var.method),
          "for variance updates\n")
      cntrl.update.var <- defaultOptimControl(var.method, lower=-5, upper=5,
                                              verbose=verbose)
    } else {
      sd.corr <- parVar$par
    }

    # overall convergence ?
    if( (parReg$rel.tol < cntrl.stop$tol) && (parVar$rel.tol < cntrl.stop$tol)
       && (parReg$convergence==0) && (parVar$convergence==0) ) 
      convergence <- 0

    # exit loop if no more change in parameters (maybe false convergence)
    if (parReg$rel.tol == 0 && parVar$rel.tol == 0)
        break
  }

  if(verbose > 0) {
    cat("\n")
    cat(as.character(Sys.time()), ":", if (convergence==0)
        "Algorithm converged" else "Algorithm DID NOT CONVERGE", "\n\n")
  }
  
  ll <- penLogLik(theta=theta,sd.corr=sd.corr,model=model)
  fisher <- penFisher(theta=theta,sd.corr=sd.corr,model=model)
  fisher.var <- marFisher(sd.corr=sd.corr, theta=theta, model=model,
                          fisher.unpen=fisher.unpen)
  
  list(fixef=head(theta,dimFE.d.O), ranef=tail(theta,dimRE), sd.corr=sd.corr,
       loglik=ll, fisher=fisher, fisherVar=fisher.var,
       convergence=convergence, dim=c(fixed=dimFE.d.O,random=dimRE))
}


#####################
# x - initial parameter values
# control arguments:
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
newtonRaphson <- function(x,fn,..., control=list(), verbose=FALSE){

  # set default values
  control.default <- list(scoreTol=1e-5, paramTol=1e-8, F.inc=0.01,
                          stepFrac=0.5, niter=30)
  control <- modifyList(control.default, control)
  
  # number of step reductions, not positive definite Fisher matrices during iterations
  steph <- notpd <- 0
  convergence <- 99
  i <- 0
  
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
                              period=52,
                              timevar="t"
                              ){
  # return formula as is if S = 0
  if(max(S) == 0) return(f)
  
  f <- paste(deparse(f), collapse="")
  # create formula
  if(length(S)==1 && S>0){
    for(i in 1:S){
      f <- paste0(f,
                  " + sin(",2*i,"*pi*",timevar,"/",period,")",
                  " + cos(",2*i,"*pi*",timevar,"/",period,")")
    }
  } else {
    for(i in 1:max(S)){
      which <- paste(i <= S,collapse=",")
      f <- paste0(f,
                  " + fe( sin(",2*i,"*pi*",timevar,"/",period,"), which=c(",which,"))",
                  " + fe( cos(",2*i,"*pi*",timevar,"/",period,"), which=c(",which,"))")
    }
  }
  return(as.formula(f, env=.GlobalEnv))
}


oneStepAhead <- function(result, # result of call to hhh4
                         tp,     # one-step-ahead predictions for time points (tp+1):nrow(stsObj)
                         verbose=TRUE,  # FIXME: this is currently unused
                         keep.estimates = FALSE){
  
  stsObj <- result$stsObj
  control <- result$control
  # get model terms
  model <- result[["terms"]]
  if(is.null(model)){
      model <- interpretControl(control,stsObj)
  }
  
  # which model: negbin or poisson?
  dimOverdisp <- model$nOverdisp
  negbin <- dimOverdisp>0

  nTime <- nrow(stsObj)
  pred <- matrix(NA_real_,nrow=length(tp:nTime)-1,ncol=ncol(stsObj))
  psi <- matrix(NA_real_,nrow=length(tp:nTime)-1,ncol=ifelse(dimOverdisp>1,ncol(stsObj),1))
    
  res <- resN <- result
  
  which <- 1
  search2 <- function(i,which){
    control$subset <- 2:i
    
    if(which==1){
      #starting values of previous time point i-1
      control$start <- list(fixed=fixef(res.old),random=ranef(res.old),sd.corr=getSdCorr(res.old))
    } else {
      #starting values of last time point
      control$start <- list(fixed=fixef(resN),random=ranef(resN),sd.corr=getSdCorr(resN))
    }
    res <- hhh4(stsObj,control)  

    return(res)
  }
  
  if(keep.estimates){
	# save value of log-likelihood (pen+mar) and parameter estimates
	params <- matrix(NA_real_,nrow=length(tp:nTime)-1,ncol=length(coef(res))+2)
	# save values of estimated covariance
	vars <- matrix(NA_real_,nrow=length(tp:nTime)-1,ncol=model$nSigma)
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

