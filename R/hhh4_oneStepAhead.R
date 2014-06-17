################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Compute one-step-ahead predictions (means) at a series of time points
###
### Copyright (C) 2011-2014 Michaela Paul and Sebastian Meyer
### $Revision: 802 $
### $Date: 2014-02-26 21:01:45 +0100 (Wed, 26 Feb 2014) $
################################################################################


oneStepAhead <- function(result, # hhh4-object (i.e. a hhh4 model fit)
                         tp,     # scalar: one-step-ahead predictions for time
                                 # points (tp+1):nrow(stsObj), or tp=c(from, to)
                         type = c("rolling", "first", "final"),
                         which.start = c("current", "final"), #if type="rolling"
                         keep.estimates = FALSE,
                         verbose = TRUE, # verbose-1 is used as verbose setting
                                         # for sequentially refitted hhh4 models
                         cores = 1) # if which.start="final", the predictions
                                    # can be computed in parallel
{
    stopifnot(inherits(result, c("ah4", "hhh4")))
    type <- match.arg(type)
    which.start <- if (type == "rolling") match.arg(which.start) else "final"
    if (cores > 1 && which.start == "current")
        stop("no parallelization for \"rolling\" if 'which.start=\"current\"'")
    startfinal <- hhh4coef2start(result)

    ## get model terms
    model <- result[["terms"]]
    if (is.null(model))
        model <- result$terms <- with(result, interpretControl(control, stsObj))
    nTime <- model$nTime
    nUnits <- model$nUnits
    psiIdx <- model$nFE + model$nd + seq_len(model$nOverdisp)
    withPsi <- length(psiIdx) > 0L
    
    ## check that tp is within the time period of the data
    maxlag <- if (is.null(result$lags) || all(is.na(result$lags)))
        1L else max(result$lags, na.rm=TRUE)
    stopifnot(tp %in% seq.int(maxlag,nTime-1L), length(tp) %in% 1:2)
    if (length(tp) == 1) tp <- c(tp, max(model$subset)-1)
    tps <- tp[1]:tp[2]
    ntps <- length(tps)
    observed <- model$response[tps+1,,drop=FALSE]
    rownames(observed) <- tps+1

    ## adjust verbosity for model refitting
    verbose <- as.integer(verbose)
    result$control$verbose <- max(0, verbose - (ntps>1))
    if (type != "rolling" && verbose > 1L) verbose <- 1L
    do_pb <- verbose == 1L
    
    ## initial fit
    fit <- if (type == "first") {
        if (do_pb)
            cat("\nRefitting model at first time point t =", tps[1L], "...\n")
        update.hhh4(result, subset.upper = tps[1L], start = startfinal,
                    keep.terms = TRUE) # need "model" -> $terms
    } else result
    if (!fit$convergence) stop("initial fit did not converge")

    ## result templates (named and filled with NA's)
    pred <- matrix(NA_real_, nrow=ntps, ncol=nUnits,
                   dimnames=list(tps+1, colnames(observed)))
    if (withPsi)
        psi <- matrix(NA_real_, nrow=ntps, ncol=length(psiIdx),
                      dimnames=list(tps, names(model$initialTheta)[psiIdx]))
    if (keep.estimates) {
        coefficients <- matrix(NA_real_,
                               nrow=ntps, ncol=length(model$initialTheta),
                               dimnames=list(tps, names(model$initialTheta)))
        Sigma.orig <- matrix(NA_real_, nrow=ntps, ncol=model$nSigma,
                             dimnames=list(tps, names(result$Sigma.orig)))
        logliks <- matrix(NA_real_, nrow=ntps, ncol=2L,
                          dimnames=list(tps, c("loglikelihood", "margll")))
    }

    ## extract predictions and stuff for specific tp from fit
    getPreds <- function (fit, tp) {
        coefs <- unname(fit$coefficients)
        c(list(pred = as.vector(
               meanHHH(coefs, fit$terms, subset=tp+1L, total.only=TRUE))),
          if (withPsi) list(psi = coefs[psiIdx]),
          if (keep.estimates) list(
              coefficients=coefs,
              Sigma.orig=unname(fit$Sigma.orig),
              logliks=c(fit$loglikelihood, fit$margll))
          )
    }

    ## compute the predictions and save
    ## pred, psi, coefficients, Sigma.orig, and logliks
    if (cores > 1L) {

        ## return value template (unnamed NA vectors)
        resTemplate <- lapply(getPreds(fit, tps[1L]), "is.na<-", TRUE)

        ## run parallel
        res <- parallel::mclapply(tps, function (tp) {
            if (verbose)
                cat("One-step-ahead prediction @ t =", tp, "...\n")
            if (type == "rolling") { # update fit
                fit <- update.hhh4(result, subset.upper=tp, start=startfinal,
                                   verbose=FALSE, # chaotic in parallel
                                   keep.terms=TRUE) # need "model" -> $terms
                if (!fit$convergence) {
                    cat("WARNING: No convergence @ t =", tp, "!\n")
                    return(resTemplate)
                }
            }
            getPreds(fit, tp)
        }, mc.preschedule=TRUE, mc.cores=cores)

        ## gather results
        .extractFromList <- function (what)
            t(vapply(res, "[[", resTemplate[[what]], what, USE.NAMES=FALSE))
        pred[] <- .extractFromList("pred")
        if (withPsi)
            psi[] <- .extractFromList("psi")
        if (keep.estimates) {
            coefficients[] <- .extractFromList("coefficients")
            Sigma.orig[] <- .extractFromList("Sigma.orig")
            logliks[] <- .extractFromList("logliks")
        }
        
    } else { ## sequential one-step ahead predictions

        if (do_pb) pb <- txtProgressBar(min=0, max=ntps, initial=0, style=3)
        for(i in seq_along(tps)) {
            if (verbose > 1L) {
                cat("\nOne-step-ahead prediction @ t =", tps[i], "...\n")
            } else if (do_pb) setTxtProgressBar(pb, i)
            
            if (type == "rolling") { # update fit
                fit.old <- fit # backup
                fit <- update.hhh4(result, subset.upper=tps[i],
                                   start=switch(which.start,
                                                current=hhh4coef2start(fit),
                                                final=startfinal),
                                   keep.terms=TRUE) # need "model" -> $terms
                if (!fit$convergence) {
                    if (do_pb) cat("\n")
                    cat("WARNING: No convergence @ t =", tps[i], "!\n")
                    ## FIXME: do a grid search ?
                    fit <- fit.old
                    next
                }
            }

            res <- getPreds(fit, tps[i])
            
            ## gather results
            pred[i,] <- res$pred
            if (withPsi)
                psi[i,] <- res$psi
            if (keep.estimates) {
                coefficients[i,] <- res$coefficients
                Sigma.orig[i,] <- res$Sigma.orig
                logliks[i,] <- res$logliks
            }
        }
        if (do_pb) close(pb)

    }

    ## done
    c(list(pred = pred, observed = observed,
           psi = if (withPsi) psi else NULL,
           allConverged = all(!is.na(pred))),
      if (keep.estimates) list(coefficients = coefficients,
                               Sigma.orig = Sigma.orig,
                               logliks = logliks)
      )
}
