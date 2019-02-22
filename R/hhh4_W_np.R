################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Non-parametric specification of neighbourhood weights in hhh4()
###
### Copyright (C) 2014,2018 Sebastian Meyer
### $Revision: 2246 $
### $Date: 2018-11-22 15:10:21 +0100 (Thu, 22. Nov 2018) $
################################################################################


### non-parametric estimation of weight function, i.e., provide each
### neighbourhood order (including 0 if from0=TRUE) up to 'maxlag' with its
### own (unconstrained) weight. For identifiability:
### - lowest order is fixed to weight=1
### - usually maxlag < max(nborder) (since only few pairs with highest orders),
###   and 'truncate' indicates if there should be zero weight for orders above
###   'maxlag' (default), or the same as for order 'maxlag'
### Thus, if from0, the parameters refer to lags 1:maxlag, otherwise 2:maxlag

W_np <- function (maxlag, truncate = TRUE, normalize = TRUE,
                  initial = log(zetaweights(2:(maxlag+from0))), from0 = FALSE,
                  to0 = truncate)  # 'to0' has been renamed to 'truncate'
{
    if (missing(maxlag)) {
        stop("'maxlag' must be specified (usually < max. neighbourhood order)")
    } else {
        stopifnot(isScalar(maxlag), maxlag >= 2 - from0) # at least one parameter
    }
    stopifnot(is.vector(initial, mode = "numeric"),
              length(initial) == maxlag + from0 - 1)

    if (!missing(to0)) {
        .Deprecated(msg = "argument 'to0' has been renamed; use 'truncate'")
        truncate <- to0
    }

    ## auxiliary expression used in 'dw' and 'd2w' below
    indicatormatrixExpr <- if (truncate) {
        quote(nbmat==nbOrder)
    } else {
        if (from0) { # maxlag = npars
            quote(if(nbOrder==npars) nbmat>=nbOrder else nbmat==nbOrder)
        } else { # maxlag = 1 + npars
            quote(if(nbOrder==1L+npars) nbmat>=nbOrder else nbmat==nbOrder)
        }
    }

    ## weights as a function of parameters and a matrix of neighbourhood orders
    w <- function (logweights, nbmat, ...) {}
    body(w) <- substitute(
    {
        weights <- exp(logweights)      # values for orders (2-from0):maxlag
        npars <- length(weights)
        W <- .WEIGHTS[1L+nbmat]         # substituted depending on 'from0'
        ## repeat last coefficient for higher orders without separate estimate
        W[is.na(W)] <- .HOWEIGHT        # substituted depending on 'truncate'
        dim(W) <- dimW <- dim(nbmat)    # nUnits x nUnits
        dimnames(W) <- dimnames(nbmat)
        .RETVAL                         # substituted depending on 'normalize'
    }, list(
        .WEIGHTS = if (from0) quote(c(1, weights)) else quote(c(0, 1, weights)),
        .HOWEIGHT = if (truncate) 0 else quote(weights[npars]),
        .RETVAL = if (normalize)
        quote(W / (norm <- .rowSums(W, dimW[1L], dimW[2L]))) else quote(W)
        ))

    ## version of w with assignment of its return value (for use in normalized
    ## versions of dw and d2w)
    .w <- w
    body(.w)[[length(body(.w))]] <-
        substitute(Wnorm <- x, list(x=body(.w)[[length(body(.w))]]))

    ## derivative of w(logweights) -> a list of matrices (one for each param.)
    if (normalize) {
        dw <- .w
        ## append code to calculate first derivatives
        body(dw) <- as.call(c(as.list(body(dw)), eval(substitute(
            expression(
                FUN <- function (nbOrder, weight) {
                    ind <- .INDICATORMATRIX
                    (ind - Wnorm*.rowSums(ind,dimW[1L],dimW[2L])) * weight/norm
                },
                mapply(FUN, .LAGS, weights, SIMPLIFY=FALSE, USE.NAMES=FALSE)
            ),
            list(.INDICATORMATRIX = indicatormatrixExpr,
                 .LAGS = if (from0) quote(seq_len(npars)) else quote(1L + seq_len(npars)))
            ))))
    } else {
        dw <- function (logweights, nbmat, ...) {}
        body(dw) <- substitute(
        {
            weights <- exp(logweights)
            npars <- length(weights)
            FUN <- function (nbOrder, weight)
                weight * (.INDICATORMATRIX)
            mapply(FUN, .LAGS, weights, SIMPLIFY=FALSE, USE.NAMES=FALSE)
        },
        list(.INDICATORMATRIX = indicatormatrixExpr,
             .LAGS = if (from0) quote(seq_len(npars)) else quote(1L + seq_len(npars))))
    }


    ## result of d2w must be a list of matrices of length npars*(npars+1L)/2L
    if (normalize) {
        d2w <- .w
        body(d2w) <- as.call(c(as.list(body(d2w)), eval(substitute(
            expression(
                seqnpars <- seq_len(npars),
                inds <- lapply(.LAGS, function (nbOrder) {
                    ind <- .INDICATORMATRIX
                    indrs <- .rowSums(ind, dimW[1L], dimW[2L])
                    list(indterm = ind - Wnorm * indrs,
                         indrs = indrs)
                }),
                k <- rep.int(seqnpars, npars), # row index
                l <- rep.int(seqnpars, rep.int(npars,npars)), # column index
                ##<- 12x faster than expand.grid(seqnpars,seqnpars)
                lowertri <- k >= l,
                ##<- and 2.5x faster than
                ##kl <- which(lower.tri(matrix(,npars,npars), diag=TRUE), arr.ind=TRUE)
                norm2 <- norm^2,
                mapply(function (k, l) weights[k] / norm2 *
                       if (k==l) {
                           inds[[k]][[1L]] * (norm - 2*weights[k]*inds[[k]][[2L]])
                       } else {
                           -weights[l] * (inds[[k]][[1L]] * inds[[l]][[2L]] +
                                          inds[[l]][[1L]] * inds[[k]][[2L]])
                       },
                       k[lowertri], l[lowertri], # inds[k[lowertri]], inds[l[lowertri]],
                       SIMPLIFY=FALSE, USE.NAMES=FALSE)
                ),
            list(.INDICATORMATRIX = indicatormatrixExpr,
                 .LAGS = if (from0) quote(seqnpars) else quote(1L + seqnpars))
            ))))
    } else { # for k=k', second derivative = first derivative, otherwise 0
        d2w <- dw
        if (length(initial) > 1) {
            ## add assignment for the return value of dw
            body(d2w)[[length(body(d2w))]] <-
                substitute(dW <- x, list(x=body(d2w)[[length(body(d2w))]]))
            ## append code to generate the list of second derivatives
            body(d2w) <- as.call(c(as.list(body(d2w)), expression(
                d2wlength <- (npars^2+npars)/2,
                ## indices of diagonal elements in x[lower.tri(x,diag=TRUE)]
                d2wdiag <- c(1L,1L+cumsum(seq.int(npars,2L))),
                d2wlist <- rep.int(list(0*nbmat), d2wlength),
                d2wlist[d2wdiag] <- dW,
                d2wlist
                )))
        }
    }

    ## Done
    environment(w) <- environment(dw) <- environment(d2w) <- .GlobalEnv
    list(w = w, dw = dw, d2w = d2w, initial = initial)
}
