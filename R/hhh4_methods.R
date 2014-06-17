################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Standard methods for hhh4-fits
###
### Copyright (C) 2010-2014 Michaela Paul and Sebastian Meyer
### $Revision: 833 $
### $Date: 2014-03-12 00:46:56 +0100 (Wed, 12 Mar 2014) $
################################################################################


## NOTE: we also apply print.hhh4 in print.summary.hhh4()
print.hhh4 <- function (x, digits = max(3, getOption("digits")-3), ...)
{
    if (!x$convergence) {
        cat('Results are not reliable! Try different starting values.\n')
        return(invisible(x))
    }
    if (!is.null(x$call)) {
        cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
            "\n\n", sep = "")
    }
    if (x$dim["random"] > 0) {
        cat('Random effects:\n')
        .printREmat(if (is.null(x$REmat)) .getREmat(x) else x$REmat,
                    digits = digits)
        cat("\nFixed effects:\n")
    } else if (x$dim["fixed"] > 0) {
        cat("Coefficients:\n")
    }
    if (x$dim["fixed"] > 0) {
        print.default(
            format(if (is.null(x$fixef)) fixef.hhh4(x, ...) else x$fixef,
                   digits=digits),
            quote = FALSE, print.gap = 2)
    } else cat("No coefficients\n")
    cat("\n")
    invisible(x)
}

## get estimated covariance matrix of random effects
.getREmat <- function (object)
{
    ## return NULL if model has no random effects
    if (is.null(REmat <- object$Sigma)) return(NULL)

    ## hhh4()$Sigma is named since r791 only -> derive names from Sigma.orig
    if (is.null(dimnames(REmat)))
        dimnames(REmat) <- rep.int(
            list(sub("^sd\\.", "",
                     names(object$Sigma.orig)[seq_len(nrow(REmat))])), 2L)
    
    attr(REmat, "correlation") <- cov2cor(REmat)
    attr(REmat, "sd") <- sqrt(diag(REmat))
    REmat
}
.printREmat <- function (REmat, digits = 4)
{
    V <- round(diag(REmat), digits=digits)
    corr <- round(attr(REmat, "correlation"), digits=digits)
    corr[upper.tri(corr,diag=TRUE)] <- ""
    V.corr <- cbind(V, corr, deparse.level=0)
    colnames(V.corr) <- c("Var", "Corr", rep.int("", ncol(corr)-1L))
    print.default(V.corr, quote=FALSE)
}

summary.hhh4 <- function (object, maxEV = FALSE, ...)
{
    ## do not summarize results in case of non-convergence
    if (!object$convergence) {
        cat('Results are not reliable! Try different starting values.\n')
	return(invisible(object))
    }
    ret <- c(object[c("call", "convergence", "dim", "loglikelihood", "margll",
                      "lags", "nTime", "nUnit")],
             list(fixef = fixef.hhh4(object, se=TRUE, ...),
                  ranef = ranef.hhh4(object, ...),
                  REmat = .getREmat(object),
                  AIC   = AIC(object),
                  BIC   = BIC(object),
                  maxEV_range = if (maxEV) unique(range(getMaxEV(object)))))
    class(ret) <- "summary.hhh4"
    return(ret)
}

print.summary.hhh4 <- function (x, digits = max(3, getOption("digits")-3), ...)
{
    ## x$convergence is always TRUE if we have a summary
    print.hhh4(x) # also works for summary.hhh4-objects

    if (!is.null(x$maxEV_range))
        cat("Epidemic dominant eigenvalue: ",
            paste(round(x$maxEV_range,2), collapse = " -- "), "\n\n")
    if(x$dim["random"]==0){
        cat('Log-likelihood:  ',round(x$loglikelihood,digits=digits-2),'\n')  
        cat('AIC:             ',round(x$AIC,digits=digits-2),'\n')
        cat('BIC:             ',round(x$BIC,digits=digits-2),'\n\n')
    } else {
        cat('Penalized log-likelihood: ',round(x$loglikelihood,digits=digits-2),'\n')  
        cat('Marginal log-likelihood:  ',round(x$margll,digits=digits-2),'\n\n')        
    }
    cat('Number of units:       ', x$nUnit, '\n')
    cat('Number of time points: ', x$nTime, '\n')
    if (!is.null(x$lags)) { # only available since surveillance 1.8-0
        if (!is.na(x$lags["ar"]) && x$lags["ar"] != 1)
            cat("Non-default autoregressive lag:  ", x$lags[["ar"]], "\n")
        if (!is.na(x$lags["ne"]) && x$lags["ne"] != 1)
            cat("Non-default neighbor-driven lag: ", x$lags[["ne"]], "\n")
    }
    cat("\n")
    invisible(x)
}

terms.hhh4 <- function (x, ...)
{
    if (is.null(x$terms))
        interpretControl(x$control,x$stsObj) else x$terms
}

nobs.hhh4 <- function (object, ...) {
    if (object$convergence) object$nObs else NA_real_
}

logLik.hhh4 <- function(object, ...)
{
    val <- if (object$convergence) object$loglikelihood else {
        warning("algorithm did not converge")
        NA_real_
    }
    attr(val, "df") <- if (object$dim["random"])
        NA_integer_ else object$dim[["fixed"]]  # use "[[" to drop the name
    attr(val, "nobs") <- nobs.hhh4(object)
    class(val) <- "logLik"
    val
}

coef.hhh4 <- function(object, se=FALSE,
                      reparamPsi=TRUE, idx2Exp=NULL, amplitudeShift=FALSE, ...)
{
    if (object$control$family == "Poisson") reparamPsi <- FALSE
    coefs <- object$coefficients
    coefnames <- names(coefs)
    idx <- getCoefIdxRenamed(coefnames, reparamPsi, idx2Exp, amplitudeShift,
                             warn=!se)
    
    ## transform and rename
    if (length(idx$Psi)) {
        coefs[idx$Psi] <- exp(-coefs[idx$Psi])  # -log(overdisp) -> overdisp
        coefnames[idx$Psi] <- names(idx$Psi)
    }
    if (length(idx$toExp)) {
        coefs[idx$toExp] <- exp(coefs[idx$toExp])
        coefnames[idx$toExp] <- names(idx$toExp)
    }
    if (length(idx$AS)) {
        coefs[idx$AS] <- sinCos2amplitudeShift(coefs[idx$AS])
        coefnames[idx$AS] <- names(idx$AS)
    }
    ## set new names
    names(coefs) <- coefnames
    
    if (se) {
        cov <- vcov.hhh4(object, reparamPsi=reparamPsi, idx2Exp=idx2Exp,
                         amplitudeShift=amplitudeShift)
        cbind("Estimates"=coefs, "Std. Error"=sqrt(diag(cov)))
    } else coefs
}

vcov.hhh4 <- function (object,
                       reparamPsi=TRUE, idx2Exp=NULL, amplitudeShift=FALSE, ...)
{
    if (object$control$family == "Poisson") reparamPsi <- FALSE
    idx <- getCoefIdxRenamed(names(object$coefficients),
                             reparamPsi, idx2Exp, amplitudeShift, warn=FALSE)
    newcoefs <- coef.hhh4(object, se=FALSE, reparamPsi=reparamPsi,
                          idx2Exp=idx2Exp, amplitudeShift=amplitudeShift)

    ## Use multivariate Delta rule => D %*% vcov %*% t(D), D: Jacobian.
    ## For idx2Exp and reparamPsi, we only transform coefficients independently,
    ## i.e. D is diagonal (with elements 'd')
    d <- rep.int(1, length(newcoefs))
    if (length(idx$Psi)) # h = exp(-psi), h' = -exp(-psi)
        d[idx$Psi] <- -newcoefs[idx$Psi]
    if (length(idx$toExp)) # h = exp(coef), h' = exp(coef)
        d[idx$toExp] <- newcoefs[idx$toExp]
    ## For the amplitude/shift-transformation, D is non-diagonal
    vcov <- if (length(idx$AS)) {
        D <- diag(d, length(d))
        D[idx$AS,idx$AS] <- jacobianAmplitudeShift(newcoefs[idx$AS])
        D %*% object$cov %*% t(D)
    } else t(t(object$cov*d)*d)  # 30 times faster than via matrix products
        
    dimnames(vcov) <- list(names(newcoefs), names(newcoefs))
    vcov
}

getCoefIdxRenamed <- function (coefnames, reparamPsi=TRUE, idx2Exp=NULL,
                               amplitudeShift=FALSE, warn=TRUE)
{
    ## indexes of overdispersion parameters
    idxPsi <- if (reparamPsi) {
        idxPsi <- grep("-log(overdisp", coefnames, fixed=TRUE)
        ## change labels from "-log(overdisp.xxx)" to "overdisp.xxx"
        names(idxPsi) <- substr(coefnames[idxPsi], start=6,
                                stop=nchar(coefnames[idxPsi])-1L)
        if (length(idxPsi) == 0L) { # backward compatibility (internal psi coef
                                    # was named "overdisp" prior to r406)
            idxPsi <- grep("^overdisp", coefnames)
            names(idxPsi) <- coefnames[idxPsi]
        }
        idxPsi
    } else NULL

    ## indexes of sine-cosine coefficients
    idxAS <- if (amplitudeShift) {
        idxAS <- sort(c(grep(".sin(", coefnames, fixed=TRUE),
                        grep(".cos(", coefnames, fixed=TRUE)))
        names(idxAS) <- sub(".sin", ".A", coefnames[idxAS], fixed=TRUE)
        names(idxAS) <- sub(".cos", ".s", names(idxAS), fixed=TRUE)
        idxAS
    } else NULL

    ## indexes of coefficients to exp()-transform
    idx2Exp <- if (length(idx2Exp)) { # index sets must be disjoint
        if (length(idxOverlap <- intersect(c(idxPsi, idxAS), idx2Exp))) {
            if (warn)
                warning("following 'idx2Exp' were ignored due to overlap: ",
                        paste(idxOverlap, collapse=", "))
            idx2Exp <- setdiff(idx2Exp, idxOverlap)
        }
        if (length(idx2Exp))
            names(idx2Exp) <- paste0("exp(", coefnames[idx2Exp], ")")
        idx2Exp
    } else NULL

    ## done
    list(Psi=idxPsi, AS=idxAS, toExp=idx2Exp)
}

fixef.hhh4 <- function (object,...)
{
    if (object$dim[1L] > 0) {
        head(coef.hhh4(object, ...), object$dim[1L])
    } else NULL
}

ranef.hhh4 <- function (object, tomatrix = FALSE, ...)
{
    if (object$dim[2L] > 0){
        ranefvec <- tail(coef.hhh4(object, ...), object$dim[2L])
    } else return(NULL)
    if (!tomatrix) return(ranefvec)

    ## transform to a nUnits x c matrix (c %in% 1:3)
    model <- terms.hhh4(object)
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

## adaption of stats::confint.default authored by the R Core Team
confint.hhh4 <- function (object, parm, level = 0.95,
                          reparamPsi=TRUE, idx2Exp=NULL, amplitudeShift=FALSE,
                          ...)
{
    cf <- coef.hhh4(object, se=TRUE, reparamPsi=reparamPsi, idx2Exp=idx2Exp,
                    amplitudeShift=amplitudeShift, ...)
    pnames <- rownames(cf)
    if (missing(parm)) 
        parm <- pnames
    else if (is.numeric(parm)) 
        parm <- pnames[parm]
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- paste(format(100*a, trim=TRUE, scientific=FALSE, digits=3), "%")
    fac <- qnorm(a)
    ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
    ses <- cf[parm,2]
    ci[] <- cf[parm,1] + ses %o% fac
    ci
}


## mean predictions for a subset of 1:nrow(object$stsObj)
predict.hhh4 <- function(object, newSubset = object$control$subset,
                         type = "response", ...)
{
    if (type == "response" &&
        all((m <- match(newSubset, object$control$subset, nomatch=0L)) > 0)) {
        ## we can extract fitted means from object
        object$fitted.values[m,,drop=FALSE]
    } else { ## means for time points not fitted (not part of object$control$subset)
        predicted <- meanHHH(coef.hhh4(object, reparamPsi=FALSE),
                             terms.hhh4(object),
                             subset=newSubset)
        if (type=="response") predicted$mean else {
            type <- match.arg(type, names(predicted))
            predicted[[type]]
        }
    }
}


### refit hhh4-model
## ...: arguments modifying the original control list
## subset.upper: refit on a subset of the data up to that time point
## use.estimates: use fitted parameters as new start values
##                (only applicable if same model)

update.hhh4 <- function (object, ..., subset.upper=NULL, use.estimates=FALSE)
{
    control <- object$control
    control <- modifyList(control, list(...))
    if (isScalar(subset.upper))
        control$subset <- control$subset[control$subset <= subset.upper]
    if (use.estimates)
        control$start <- hhh4coef2start(object)
    hhh4(object$stsObj, control)
}


## convert fitted parameters to a list suitable for control$start
hhh4coef2start <- function (fit)
    list(fixed = fit$coefficients[seq_len(fit$dim[1L])],
         random = fit$coefficients[fit$dim[1L]+seq_len(fit$dim[2L])],
         sd.corr = fit$Sigma.orig)

## character vector of model components that are "inModel"
componentsHHH4 <- function (object)
    names(which(sapply(object$control[c("ar", "ne", "end")], "[[", "inModel")))

## deviance residuals
residuals.hhh4 <- function (object, type = c("deviance", "response"), ...)
{
    type <- match.arg(type)
    obs <- observed(object$stsObj)[object$control$subset,]
    fit <- fitted(object)
    if (type == "response")
        return(obs - fit)

    ## deviance residuals
    ## Cf. residuals.ah, it calculates:
    ## deviance = sign(y - mean) * sqrt(2 * (distr(y) - distr(mean)))
    ## pearson = (y - mean)/sqrt(variance)
    dev.resids <- switch(object$control$family,
        "Poisson" = poisson()$dev.resids,
        "NegBin1" = {
            psi <- exp(object$coefficients[object$dim[1L]])
            negative.binomial(psi)$dev.resids
        },
        "NegBinM" = {
            psicoefidx <- seq.int(to=object$dim[1L],
                                  length.out=object$nUnit)
            psi <- matrix(
                exp(object$coefficients[psicoefidx]),
                object$nTime, object$nUnit, byrow=TRUE)
            negative.binomial(psi)$dev.resids # CAVE: non-standard use
        },
        stop("not implemeted for \"", object$control$family, "\"-models"))

    di2 <- dev.resids(y=obs, mu=fit, wt=1)
    sign(obs-fit) * sqrt(pmax.int(di2, 0))
}
