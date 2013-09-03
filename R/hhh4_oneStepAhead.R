################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Compute one-step-ahead predictions (means) at a series of time points
###
### Copyright (C) 2011-2013 Michaela Paul and Sebastian Meyer
### $Revision: 594 $
### $Date: 2013-07-08 16:05:55 +0200 (Mon, 08 Jul 2013) $
################################################################################


oneStepAhead <- function(result, # ah4-object (i.e. a hhh4 model fit)
                         tp,     # scalar: one-step-ahead predictions for time
                                 # points (tp+1):nrow(stsObj), or tp=c(from, to)
                         which.start = c("current", "final", "none"),
                         keep.estimates = FALSE,
                         verbose = TRUE) # verbose-1 is used as verbose setting
                                         # for sequentially refitted hhh4 models
{
    stopifnot(inherits(result, "ah4"))
    which.start <- match.arg(which.start)
    refit <- which.start != "none"
    use.current <- which.start == "current"
    ## i.e., use fitted parameters from previous time point as initial values
    startfinal <- ah4coef2start(result)
    
    ## get model terms
    model <- result[["terms"]]
    if (is.null(model))
        model <- result$terms <- with(result, interpretControl(control, stsObj))
    nTime <- model$nTime
    nUnits <- model$nUnits
    
    ## check that tp is within the time period of the data
    stopifnot(tp %in% seq_len(nTime-1L), length(tp) %in% 1:2)
    if (length(tp) == 1) tp <- c(tp, nTime-1) # historical default
    tps <- tp[1]:tp[2]
    ntps <- length(tps)
    observed <- model$response[tps+1,,drop=FALSE]
    rownames(observed) <- tps+1
    
    ## adjust verbosity for model refitting
    result$control$verbose <- verbose - (ntps>1)

    ## initialize result
    pred <- matrix(NA_real_, nrow=ntps, ncol=nUnits,
                   dimnames=list(tps, colnames(observed)))
    psi <- if (model$nOverdisp > 0) {
        psiNames <- grep("overdisp", names(model$initialTheta), value=TRUE)
        matrix(NA_real_, nrow=ntps, ncol=model$nOverdisp,
               dimnames=list(tps, psiNames))
    } else NULL
    if (keep.estimates) {
	coefficients <- matrix(NA_real_,
                               nrow=ntps, ncol=length(model$initialTheta),
                               dimnames=list(tps, names(model$initialTheta)))
	Sigma.orig <- matrix(NA_real_, nrow=ntps, ncol=model$nSigma,
                             dimnames=list(tps, names(result$Sigma.orig)))
        logliks <- matrix(NA_real_, nrow=ntps, ncol=2L,
                          dimnames=list(tps, c("loglikelihood", "margll")))
    } else {
        coefficients <- Sigma.orig <- logliks <- NULL
    }

    ## sequential one-step ahead predictions
    fit <- result
    for(i in seq_along(tps)) {
        if (verbose) cat(tps[i], "\n")
        if (refit) {
            fit.old <- fit
            fit <- update.ah4(result, subset.upper=tps[i],
                              start=if (use.current)
                                    ah4coef2start(fit.old) else startfinal,
                              keep.terms=TRUE) # need "model" -> $terms
        }
        if (fit$convergence) {
            coefs <- coef(fit, reparamPsi=FALSE)
            pred[i,] <- meanHHH(coefs, fit$terms,
                                subset=tps[i]+1, total.only=TRUE)
            if (model$nOverdisp > 0)
                psi[i,] <- coefs[psiNames]
            if (keep.estimates) {
                coefficients[i,] <- coefs
                Sigma.orig[i,] <- getSdCorr(fit)
                logliks[i,] <- c(fit$loglikelihood, fit$margll)
            }
        } else { # do a grid search ?
            if (verbose) cat("-> NO convergence at time point", tps[i], "\n")
            if (refit) fit <- fit.old
        }
    }
    
    list(pred=pred, psi=psi, observed=observed, allConverged=all(!is.na(pred)),
         coefficients=coefficients, Sigma.orig=Sigma.orig, logliks=logliks)
}
