print.ah4 <- function (x, digits = max(3, getOption("digits") - 3),reparamPsi=TRUE, ...) 
{
  if(!x$convergence)
    cat('Results are not reliable! Try different starting values. \n')
  else {
    if(!is.null(x$call)){
      cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
	}

    if (length(ranef(x))) {
      cat('Random effects: \n')
      V <- as.matrix(round(diag(getCov(x)),digits=digits))
      corr <- round(getCov(x),digits=digits)
      corr[upper.tri(corr,diag=TRUE)] <- ""
      V.corr <- cbind(V,corr)
      colnames(V.corr) <- c("Var","Corr",rep("",ncol(V.corr)-2))
      print(noquote(V.corr))

	  cat("\nFixed effects:\n")
	  print.default(format(fixef(x, reparamPsi=TRUE,...), digits = digits), print.gap = 2, 
		  quote = FALSE)
    } else if (length(fixef(x))) {
        cat("Coefficients:\n")
        print.default(format(fixef(x, reparamPsi=TRUE,...), digits = digits), print.gap = 2, 
            quote = FALSE)
    }
    else cat("No coefficients\n")
    cat("\n")
    invisible(x)
  }
}

summary.ah4 <- function(object, ...){
  # do not summarize results in case of non-convergence
  if(!object$convergence){
    cat('Results are not reliable! Try different starting values. \n')
	return(invisible(object))
  }

  ret <- object[c("call","convergence","dim","loglikelihood","margll","nTime","nUnit")]

  # get estimated covariance matrix of random effects
  if(object$dim["random"]>0){
	REmat <- object$Sigma
	attr(REmat, "correlation") <- cov2cor(REmat)
	attr(REmat, "sd") <- sqrt(diag(REmat))
    rownames(REmat) <- colnames(REmat) <- gsub("sd.", "", names(object$Sigma.orig))[grep("sd.", names(object$Sigma.orig))]

	aic <- NULL
	bic <- NULL
  } else{
	REmat <- NULL
	aic <- AIC(object)
    bic <- AIC(object,k=log(object$nObs))
  }

  fe <- fixef(object,se=TRUE, ...)
  re <- ranef(object)

  ret <- c(ret, list(fixef=fe, ranef=re, REmat=REmat, AIC=aic, BIC=bic))
  class(ret) <- "summary.ah4"
  return(ret)
}

print.summary.ah4 <- function(x, digits = max(3, getOption("digits") - 3), ...){
  if(!x$convergence){
    cat('Results are not reliable! Try different starting values. \n')
	invisible(x)
  } else {
    if(!is.null(x$call)){
      cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    }
    
    if(x$dim["random"]>0){
      cat('Random effects: \n')
      V <- as.matrix(round(diag(x$REmat),digits=digits))
      corr <- round(attr(x$REmat, "correlation"),digits=digits)
      corr[upper.tri(corr,diag=TRUE)] <- ""
      V.corr <- cbind(V,corr)
      colnames(V.corr) <- c("Var","Corr",rep("",ncol(V.corr)-2))
      print(noquote(V.corr))
    }
        
    cat('\nFixed effects: \n')
    print(round(x$fixef,digits=digits),print.gap=2)

    if(x$dim["random"]==0){
      cat('\nlog-likelihood:   ',round(x$loglikelihood,digits=digits-2),'\n')  
      cat('AIC:              ',round(x$AIC,digits=digits-2),'\n')
      cat('BIC:              ',round(x$BIC,digits=digits-2),'\n\n')
    } else {
      cat('\npenalized log-likelihood: ',round(x$loglikelihood,digits=digits-2),'\n')  
      cat('marginal log-likelihood:  ',round(x$margll,digits=digits-2),'\n\n')        
    }
    
    cat('Number of units:         ',x$nUnit,'\n')
    cat('Number of time points:   ',x$nTime,'\n\n')
    
    if(any(x$fixef[,"Estimates"]== -20)){
      cat("Warning: coefficients \'",rownames(x$fixef)[x$fixef[,"Estimates"]==-20] ,"\' could not be estimated accurately\n\n")
    }
  }
  
  invisible(x)
}

logLik.ah4 <- function(object,...){
  if(!object$convergence)
	stop("algorithm did not converge\n")

  val <- object$loglikelihood
  if(object$dim["random"]==0){
	attr(val, "df") <- length(coef(object))
	attr(val, "nobs") <- object$nObs
	class(val) <- "logLik"
  }
   return(val)
}

AIC.ah4 <- function (object, ..., k = 2)
{
    if(object$dim["random"]==0){
        NextMethod("AIC")
    } else {
	warning("AIC not well defined for models with random effects")
        NA_real_
    }
}

coef.ah4 <- function(object, se=FALSE, reparamPsi=TRUE, idx2Exp=NULL,
                     amplitudeShift=FALSE, ...)
{
  coefs <- object$coefficients
  stdErr <- object$se

  if(reparamPsi && object$control$family!="Poisson"){
    #extract psi coefficients
    index <- grep("-log(overdisp",names(coefs), fixed=TRUE)

    if (length(index) == 0L) { # backward compatibility (internal psi coef
                               # was named "overdisp" prior to r406)
      index <- grep("^overdisp", names(coefs))
    } else {
      psi.names <- names(coefs)[index]
      # change labels from "-log(overdisp.xxx)" to "overdisp.xxx"
      names(coefs)[index] <- substr(psi.names, start=6, stop=nchar(psi.names)-1)
    }
    
    #transform psi coefficients
    coefs[index] <- exp(-coefs[index])
    # se's using Delta rule: se[h(psi)] = se[psi] * |h'(psi)|
    # h = exp(-psi), h' = -exp(-psi)
    D <- diag(coefs[index],length(index))
    stdErr[index] <- sqrt(diag(D %*% object$cov[index,index] %*% t(D)))
  }
  
  if(!is.null(idx2Exp)){
    # extract coefficients on log-scale
    exp.names <- names(coefs)[idx2Exp]
    # change labels
    names(coefs)[idx2Exp] <- paste("exp(",exp.names,")",sep="")
    
    # transform
    coefs[idx2Exp] <- exp(coefs[idx2Exp])
    D <- diag(coefs[idx2Exp],length(idx2Exp))
    stdErr[idx2Exp] <- sqrt(diag(D %*% object$cov[idx2Exp,idx2Exp] %*% t(D)))
  }
  
  
  if(amplitudeShift){
    indexAS <- sort(c(grep(c(".sin"),names(coefs),fixed=TRUE),grep(c(".cos"),names(coefs),fixed=TRUE)))
    namesSinCos <- names(coefs)[indexAS]
    namesSinCos <- gsub(".sin",".A",namesSinCos)
    namesSinCos <- gsub(".cos",".s",namesSinCos)
    names(coefs)[indexAS] <- namesSinCos
    coefs[indexAS] <- sinCos2amplitudeShift(coefs[indexAS])
    D <- jacobianAmplitudeShift(coefs[indexAS])
    cov <- D %*% object$cov[indexAS,indexAS] %*% t(D)
    stdErr[indexAS] <- sqrt(diag(cov))
  }

  if(se)
    return(cbind("Estimates"=coefs,"Std. Error"=stdErr))
  else
    return(coefs)
}

fixef.ah4 <- function(object,...){
  if(object$dim[1]>0){
    return(head(coef(object,...), object$dim[1]))
  } else return(NULL)
}

ranef.ah4 <- function(object, tomatrix = FALSE, ...){
  if(object$dim[2]>0){
    ranefvec <- tail(coef(object,...), object$dim[2])
  } else return(NULL)
  if (!tomatrix) return(ranefvec)

  ## transform to a nUnits x c matrix (c %in% 1:3)
  model <- interpretControl(object$control,object$stsObj)
  idxRE <- model$indexRE
  idxs <- unique(idxRE)
  names(idxs) <- model$namesFE[idxs]
  mat <- sapply(idxs, function (idx) {
      RE <- ranefvec[idxRE==idx]
      Z <- model$terms["Z.intercept",][[idx]]
      "%m%" <- get(model$terms["mult",][[idx]])
      Z %m% RE
  })
  rownames(mat) <- colnames(model$response)
  return(mat)
}

confint.ah4 <- function (object, parm, level = 0.95, reparamPsi = TRUE, idx2Exp = NULL, amplitudeShift = FALSE, ...) 
{
    cf <- coef(object, reparamPsi=reparamPsi, idx2Exp=idx2Exp, amplitudeShift=amplitudeShift, se=TRUE)
    pnames <- rownames(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100 * a, trim = TRUE, scientific = FALSE, digits = 3),"%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ses <- cf[parm,2]
    ci[] <- cf[parm,1] + ses %o% fac
    ci
}

## mean predictions for a subset of 1:nrow(object$stsObj)
predict.ah4 <- function(object, newSubset=control$subset, type="response", ...)
{
    control <- object$control
    if (type == "response" &&
        all((m <- match(newSubset, control$subset, nomatch=0L)) > 0)) {
        ## we can extract fitted means from object
        object$fitted.values[m,,drop=FALSE]
    } else { ## means for time points not fitted (not part of control$subset)
        data <- object$stsObj
        coefs <- coef(object, reparamPsi=FALSE)
        model <- interpretControl(control, data)
        predicted <- meanHHH(coefs, model, subset=newSubset)
        if (type=="response") predicted$mean else {
            type <- match.arg(type, names(predicted))
            predicted[[type]]
        }
    }
}

# plot estimated mean 
# this needs to be made more customizable
plot.ah4 <- function(x, i=1, m=NULL, ylim=NULL,
                     ylab="No. infected", xlab="", title=NULL,
                     col=c("grey30","grey60","grey85"), border=col,
                     cex=.6, pch=19, hide0s=FALSE, legend=FALSE, ...)
{
  stopifnot(length(i) == 1)
  if(is.character(i) && is.na(i <- match(.i <- i, colnames(x$stsObj))))
      stop("region '", .i, "' does not exist")
  if(is.null(title)) title <- colnames(x$stsObj)[i]
  obs <- observed(x$stsObj)[x$control$subset,i]
  if(is.null(ylim)) ylim <- c(0, max(obs,na.rm=TRUE))

  start <- x$stsObj@start
  start[2] <-  start[2]+1
  plot(ts(obs,frequency=x$stsObj@freq,start=start),
       ylim=ylim,ylab=ylab,type="n",las=1,xlab=xlab)
  title(main=title,line=0.5)

  if(is.null(m))
    m <- meanHHH(coef(x,reparamPsi=FALSE),interpretControl(x$control,x$stsObj))
  m$endemic[is.na(m$mean)] <- 0
  m$epi.own[is.na(m$mean)] <- 0
  m$mean[is.na(m$mean)] <- 0
    
  tp <- (1:(length(obs)))/x$stsObj@freq +x$stsObj@start[1] + x$stsObj@start[2]/x$stsObj@freq
  if(!is.null(dim(as.matrix(m$epidemic)))) { # FIXME: will this ever be FALSE?
    polygon(c(tp[1],tp,tail(tp,1)),c(0,as.matrix(m$mean)[,i],0),col=col[1],border=border[1])
    if(!is.null(dim(m$epi.own)))
      polygon(c(tp[1],tp,tail(tp,1)),c(0,as.matrix(m$endemic)[,i]+as.matrix(m$epi.own)[,i],0),col=col[2],border=border[2])
  }
  polygon(c(tp[1],tp,tail(tp,1)),c(0,as.matrix(m$endemic)[,i],0),col=col[3],border=border[3])
  ptidx <- if (hide0s) obs > 0 else TRUE
  points(tp[ptidx],obs[ptidx],pch=pch,cex=cex)
 
  if(legend){
      legend("topright", bty="n", legend=c("observed","ne","ar","end"),
             pch=c(pch,rep(NA,3)), lty=c(NA,rep(1,3)), lwd=8, pt.lwd=1,
             col=c(1,col))
  }

}


### refit hhh4-model
## ...: arguments modifying the original control list
## subset.upper: refit on a subset of the data up to that time point
## use.estimates: use fitted parameters as new start values
##                (only applicable if same model)

update.ah4 <- function (object, ..., subset.upper=NULL, use.estimates=FALSE)
{
    control <- object$control
    control <- modifyList(control, list(...))
    if (isScalar(subset.upper))
        control$subset <- control$subset[control$subset <= subset.upper]
    if (use.estimates)
        control$start <- ah4coef2start(object)
    hhh4(object$stsObj, control)
}


## convert fitted parameters from "ah4" to list suitable for control$start
ah4coef2start <- function (fit)
    list(fixed = fixef(fit, reparamPsi=FALSE),
         random = ranef(fit, reparamPsi=FALSE),
         sd.corr = getSdCorr(fit))
