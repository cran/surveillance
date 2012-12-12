################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Helper functions for neighbourhood weight matrices in hhh4()
###
### Copyright (C) 2012 Sebastian Meyer
### $Revision: 444 $
### $Date: 2012-10-19 15:51:41 +0200 (Fr, 19. Okt 2012) $
################################################################################


### check ne$weights specification

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

checkWeightsArray <- function (W, nUnits, nTime, name)
{
    if (!is.array(W))
        stop("'", name, "' must return a matrix/array")
    if (any(dim(W)[1:2] != nUnits) || isTRUE(dim(W)[3] != nTime))
        stop("'", name, "' must conform to dimensions ",
             nUnits, " x ", nUnits, " (x ", nTime, ")")
    if (any(!diag(W) %in% 0))
        stop("'diag(", name, ")' must only contain zeroes")
    if (any(is.na(W)))
        stop("missing values in '", name, "' are not allowed")
}

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
        if (dim.d == 1L) {
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



### calculate the weighted sum of counts of adjacent (or all other) regions
### i.e. the nTime x nUnit matrix with elements ne_ti = sum_j w_jit * y_jt
## W is either a nUnits x nUnits matrix of time-constant weights w_ji
## or a nUnits x nUnits x nTime array of time-varying weights

weightedSumNE <- function (observed, weights, lag)
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  tY <- t(observed)
  
  timeconstantweights <- length(dim(weights)) == 2L
  selecti <- if (timeconstantweights) quote(weights[,i]) else quote(weights[,i,])
    
  res <- matrix(NA_real_, nrow = nTime, ncol = nUnits,
                dimnames = list(NULL, colnames(observed)))
  for(i in seq_len(nUnits)){
    weights.i <- eval(selecti)
    weightedObs <- tY * weights.i
    res[,i] <- colSums(weightedObs, na.rm=TRUE)
  }

  reslagged <- rbind(matrix(NA_real_, lag, nUnits),
                     res[seq_len(nTime-lag),,drop=FALSE])
  reslagged
}

## slower alternative, where the weights are always converted to a 3D array
weightedSumNE.old <- function(observed, weights, lag)
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  
  nhood <- array(weights, c(nUnits,nUnits,nTime))
  
  res <- matrix(NA_real_, nrow = nTime, ncol = nUnits,
                dimnames = list(NULL, colnames(observed)))
  for(i in seq_len(nUnits)){
    weights.i <- t(nhood[,i,])
    weightedObs <- observed * weights.i
    res[,i] <- rowSums(weightedObs, na.rm=TRUE)
  }
  
  rbind(matrix(NA_real_, lag, nUnits), head(res, nTime-lag))
}






###############################################
### predefined parametric weight structures ###
###############################################


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


### Construct weight matrix wji according to the Zeta-distribution with respect
### to the orders of neighbourhood (in nbmat, as e.g. obtained from nbOrder()),
### optionally fulfilling rowSums(wji) = 1
## As a formula (for j != i, otherwise wji = 0):
## - for shared=TRUE: wji = pzeta(oji; rho, maxlag) / sum_k I(ojk == oji)
## - for shared=FALSE: wji = pzeta(oji; rho, maxlag) / sum_k pzeta(ojk; rho, maxlag)
## Here, oji = oij is the order of nb of i and j,
## and pzeta(o; rho, m) = o^-rho / sum_{r=1}^m r^-rho is the Zeta-distribution
## on 1:m (also called Zipf's law).
## For shared=TRUE and normalize=FALSE, maxlag should not be greater than
## min_j(max_i(oji)), such that every region has neighbours up to order 'maxlag'
## and higher orders can not be infected. Otherwise, regions with not as
## high-order neighbours would not sum their weights to 1 (but lower).
## For shared=FALSE, maxlag=Inf yields the weights
## wji = oji^-\rho / sum_k ojk^-\rho
## In both cases, maxlag=1 yields the classical weights wji=1/nj.

zetaweights <- function (nbmat, maxlag = Inf, rho = 1,
                         normalize = TRUE, shared = FALSE)
{
    ## check maxlag
    if (!is.finite(maxlag)) maxlag <- max(nbmat)

    ## raw (non-normalized) zeta-distribution on 1:maxlag
    zetaweights <- c(0, seq_len(maxlag)^-rho)
    ##<- first 0 is for lag 0 (for instance, diag(nbmat))

    ## replace order by zetaweight of that order
    wji <- zetaweights[nbmat + 1L]
    wji[is.na(wji)] <- 0
    dim(wji) <- dim(nbmat)
    dimnames(wji) <- dimnames(nbmat)

    if (shared) {
        ## multiplicity of orders by row: dim(nbmat)==dim(multbyrow)
        multbyrow <- t(apply(nbmat, 1, function(x) table(x)[as.character(x)]))
        ## neighbours of same order share the zetaweight for that order
        wji <- wji / sum(zetaweights) / multbyrow
    }
    if (normalize) { # normalize such that each row sums to 1
        wji <- wji / rowSums(wji)
    }

    ## Done
    wji
}



### powerlaw weights
## in the non-truncated case, i.e. maxlag = max(nbmat),
## the raw powerlaw weights are defined as w_ji = o_ji^-d,
## and with (row-)normalization we have    w_ji = o_ji^-d / sum_k o_jk^-d

powerlaw <- function (maxlag, normalize = TRUE, initial = 1)
                                        # atm, only shared=FALSE is supported
{
    if (missing(maxlag)) {
        stop("'maxlag' must be specified (e.g., maximum neighbourhood order)")
    }

    ## main function which returns the weight matrix
    weights.call <- call("zetaweights", quote(nbmat), maxlag, quote(d),
                         normalize, FALSE)
    weights <- as.function(c(alist(d=, nbmat=, ...=), weights.call),
                           envir=.GlobalEnv)

    ## construct derivatives with respect to "d"
    dweights <- d2weights <- as.function(c(alist(d=, nbmat=, ...=),
        substitute({
            W <- weights.call
            is.na(nbmat) <- nbmat == 0L # w_jj(d) = 0 => w_jj'(d) = 0, sum over j!=i
            tmp1a <- log(nbmat)
            norm <- rowSums(nbmat^-d, na.rm=TRUE) # unused for raw weights
            tmp1b <- rowSums(nbmat^-d * -log(nbmat), na.rm=TRUE)/norm # set to 0 for raw
            tmp1 <- tmp1a + tmp1b
            tmp2 <- rowSums(nbmat^-d * log(nbmat)^2, na.rm=TRUE)/norm - tmp1b^2 # for 2nd deriv
            deriv <- W * -tmp1          # for d2weights: W * (tmp1^2 - tmp2)
            deriv[is.na(deriv)] <- 0
            deriv
        }, list(weights.call=weights.call))), envir=.GlobalEnv)

    ## adaptions for dweights and d2weights
    body(dweights)[[grep("^tmp2 <-", body(dweights))]] <- NULL
    body(d2weights)[[grep("^deriv <-", body(d2weights))]] <-
        quote(deriv <- W * (tmp1^2 - tmp2))
    
    ## simplifications for raw weights
    if (!normalize) {
        body(dweights)[[grep("^norm", body(dweights))]] <-
            body(d2weights)[[grep("^norm", body(d2weights))]] <- NULL
        body(dweights)[[grep("^tmp1b <-", body(dweights))]] <-
            body(d2weights)[[grep("^tmp1b <-", body(d2weights))]] <- quote(tmp1b <- 0)
        body(d2weights)[[grep("^tmp2 <-", body(d2weights))]] <- quote(tmp2 <- 0)
    }
    
    ## return list of functions
    list(w=weights, dw=dweights, d2w=d2weights, initial=initial)
}


