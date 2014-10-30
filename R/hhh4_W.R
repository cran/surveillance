################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Helper functions for neighbourhood weight matrices in hhh4()
###
### Copyright (C) 2012-2014 Sebastian Meyer
### $Revision: 1025 $
### $Date: 2014-09-22 17:58:48 +0200 (Mon, 22 Sep 2014) $
################################################################################


checkNeighbourhood <- function (neighbourhood)
{
    ## setValidity() in sts.R only guarantees correct 'dim' and 'dimnames'
    ## we also assert numeric or logical matrix with non-NA entries
    ## FIXME: However, we currently don't check for symmetry and for zeros on
    ## the diagonal...
    stopifnot(is.matrix(neighbourhood),
              nrow(neighbourhood) == ncol(neighbourhood),
              is.numeric(neighbourhood) | is.logical(neighbourhood),
              is.finite(neighbourhood))
    invisible(TRUE)
}


### calculate the weighted sum of counts of adjacent (or all other) regions
### i.e. the nTime x nUnit matrix with elements ne_ti = sum_j w_jit * y_jt
## W is either a nUnits x nUnits matrix of time-constant weights w_ji
## or a nUnits x nUnits x nTime array of time-varying weights

weightedSumNE <- function (observed, weights, lag)
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  tY <- t(observed)                     # -> nUnits x nTime
  
  res <- apply(weights, 2L, function (wi)
               ## if dim(weights)==2 (time-constant weights), length(wi)=nUnits,
               ## if dim(weights)==3, wi is a matrix of size nUnits x nTime
               .colSums(tY * wi, nUnits, nTime, na.rm=TRUE))
  
  rbind(matrix(NA_real_, lag, nUnits),
        res[seq_len(nTime-lag),,drop=FALSE])
}



##################################
### check ne$weights specification
##################################


### checks for a fixed matrix/array

checkWeightsArray <- function (W, nUnits, nTime, name)
{
    if (!is.array(W))
        stop("'", name, "' must return a matrix/array")
    if (any(dim(W)[1:2] != nUnits) || isTRUE(dim(W)[3] != nTime))
        stop("'", name, "' must conform to dimensions ",
             nUnits, " x ", nUnits, " (x ", nTime, ")")
    if (any(is.na(W)))
        stop("missing values in '", name, "' are not allowed")
    diags <- if (is.matrix(W)) diag(W) else apply(W, 3, diag)
    if (any(diags != 0))
        warning("'", name, "' has nonzeros on the diagonal",
                if (!is.matrix(W)) "s")
}


### check parametric weights specification consisting of a list of:
## - three functions: w, dw, and d2w
## - a vector of initial parameter values

checkWeightsFUN <- function (object)
{
    fnames <- paste0(c("","d","d2"), "w")
    if (any(!sapply(object[fnames], is.function)))
        stop("parametric weights require functions ",
             paste0("'", fnames, "'", collapse=", "))
    if (any(!sapply(object[fnames], function(FUN) length(formals(FUN)) >= 3L)))
        stop("parametric weights functions must accept (not necessarily use)",
             "\n  at least 3 arguments (parameter vector, ",
             "neighbourhood order matrix, data)")
    if (!is.vector(object$initial, mode="numeric") ||
        length(object$initial) == 0L)
        stop("parametric weights require initial parameter values")
    TRUE
}


### entry function for checks in hhh4()

checkWeights <- function (weights, nUnits, nTime,
                          nbmat, data)  # only used for parametric weights
{
    name <- "control$ne$weights"

    ## check specification
    testweights <- if (is.array(weights)) weights else {
        if (is.list(weights) && checkWeightsFUN(weights)
            && checkNeighbourhood(nbmat)) {
            if (all(nbmat %in% 0:1))
                warning("'", deparse(substitute(nbmat)),
                        "' is binary (should contain",
                        " general neighbourhood orders)")
            weights$w(weights$initial, nbmat, data)
        } else {
            stop("'", name, "' must be a matrix/array or a list of functions")
        }
    }
    
    ## apply matrix/array checks
    if (is.list(weights)) { # parametric weights
        checkWeightsArray(testweights, nUnits, nTime, name=paste0(name, "$w"))
        dim.d <- length(weights$initial)
        dw <- weights$dw(weights$initial, nbmat, data)
        d2w <- weights$d2w(weights$initial, nbmat, data)
        if (dim.d == 1L && !is.list(dw) && !is.list(d2w)) {
            checkWeightsArray(dw, nUnits, nTime, name=paste0(name, "$dw"))
            checkWeightsArray(d2w, nUnits, nTime, name=paste0(name, "$d2w"))
        } else {
            if (!is.list(dw) || length(dw) != dim.d)
                stop("'", name, "$dw' must return a list (of matrices/arrays)",
                     " of length ", dim.d)
            if (!is.list(d2w) || length(d2w) != dim.d*(dim.d+1)/2)
                stop("'", name, "$d2w' must return a list (of matrices/arrays)",
                     " of length ", dim.d*(dim.d+1)/2)
            lapply(dw, checkWeightsArray, nUnits, nTime,
                   name=paste0(name, "$dw[[i]]"))
            lapply(d2w, checkWeightsArray, nUnits, nTime,
                   name=paste0(name, "$d2w[[i]]"))
        }
    } else checkWeightsArray(testweights, nUnits, nTime, name=name)
    
    ## Done
    invisible(TRUE)
}



#############################################
### Utility functions for fitted hhh4-objects
#############################################


### extract the (final) weight matrix/array from a fitted hhh4 object

getNEweights <- function (object, pars = coefW(object))
{
    neweights <- object$control$ne$weights
    
    if (!is.list(neweights))  # NULL or fixed weight structure
        return(neweights)

    ## parametric weights
    nd <- length(neweights$initial)
    if (length(pars) != nd) stop("'pars' must be of length ", nd)
    neweights$w(pars, neighbourhood(object$stsObj), object$control$data)
}


### extract parameters of neighbourhood weights from hhh4-object or coef vector

coefW <- function (object)
{
    coefs <- if (inherits(object, c("hhh4","ah4"))) object$coefficients else object
    coefs[grep("^neweights", names(coefs))]
}
