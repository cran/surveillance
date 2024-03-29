################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Scoring rules as discussed in:
### Predictive model assessment for count data
### Czado, C., Gneiting, T. & Held, L. (2009)
### Biometrics 65:1254-1261
###
### Copyright (C) 2010-2012 Michaela Paul
### Copyright (C) 2014-2015, 2017-2019, 2022 Sebastian Meyer
################################################################################


## logarithmic score
## logs(P,x) = -log(P(X=x))

.logs <- function (px)
    -log(px)

logs <- function (x, mu, size=NULL) {
    if (is.null(size)) {
        - dpois(x, lambda=mu, log=TRUE)
    } else {
        - dnbinom(x, mu=mu, size=size, log=TRUE)
    }
}


## squared error score
## ses(P,x) = (x-mu_p)^2

ses <- function (x, mu, size=NULL) {
    (x-mu)^2
}


## normalized squared error score (IMPROPER)
## nses(P,x) = ((x-mu_p)/sigma_p)^2

nses <- function (x, mu, size=NULL) {
    sigma2 <- if (is.null(size)) mu else mu * (1 + mu/size)
    ((x-mu)^2) / sigma2
}


## Dawid-Sebastiani score
## dss(P,x) = ((x-mu_p)/sigma_p)^2 + 2*log(sigma_p)

.dss <- function (meanP, varP, x)
    (x-meanP)^2 / varP + log(varP)

dss <- function (x, mu, size=NULL)
    .dss(meanP = mu,
         varP = if (is.null(size)) mu else mu * (1 + mu/size),
         x = x)


## ranked probability score
## rps(P,x) = sum_0^Kmax {P(X<=k) - 1(x<=k)}^2

## for a single prediction (general formulation)
.rps <- function (P, ..., x, kmax, kmin = 0, tolerance = sqrt(.Machine$double.eps))
{
    ## compute P(X<=k)
    k <- kmin:kmax
    Pk <- P(k, ...)

    ## check precision
    if ((1 - Pk[length(Pk)])^2 > tolerance)
        warning("finite sum approximation error larger than tolerance=",
                format(tolerance))

    ## compute the RPS
    sum((Pk - (x <= k))^2)
}

## for a single Poisson prediction
rps_1P <- function (x, mu, k=40, tolerance=sqrt(.Machine$double.eps)) {
    ## return NA for non-convergent fits (where mu=NA)
    if (is.na(x) || is.na(mu)) return(NA_real_)
    ## finite sum truncation at mean + k*sd, but not before x
    kmax <- max(x, ceiling(mu + k*sqrt(mu)))
    ## skip low k for large mu
    kmin <- if (mu > k^2) floor(mu - k*sqrt(mu)) else 0
    ## compute the RPS
    .rps(P = ppois, lambda = mu, x = x,
         kmax = kmax, kmin = kmin, tolerance = tolerance)
}

## for a single NegBin prediction
rps_1NB <- function (x, mu, size, k=40, tolerance=sqrt(.Machine$double.eps)) {
    ## return NA for non-convergent fits (where mu=NA)
    if (anyNA(c(x, mu, size))) return(NA_real_)
    ## finite sum truncation at mean + k*sd, but not before x
    sigma2 <- mu * (1 + mu/size)
    kmax <- max(x, ceiling(mu + k*sqrt(sigma2)))
    ## skip low k for large mu
    kmin <- max(0, floor(mu - k*sqrt(sigma2)))
    ## protect against wide NegBin (excessive memory consumption)
    if (kmax - kmin > 1e8) {
        warning("quasi-continuous NegBin distribution (too wide); returning NA")
        return(NA_real_)
    }
    ## compute the RPS
    .rps(P = pnbinom, mu = mu, size = size, x = x,
         kmax = kmax, kmin = kmin, tolerance = tolerance)
}

## vectorized version
rps <- function (x, mu, size=NULL, k=40, tolerance=sqrt(.Machine$double.eps)) {
    res <- if (is.null(size)) {
        mapply(rps_1P, x=x, mu=mu,
               MoreArgs=list(k=k, tolerance=tolerance),
               SIMPLIFY=TRUE, USE.NAMES=FALSE)
    } else {
        mapply(rps_1NB, x=x, mu=mu, size=size,
               MoreArgs=list(k=k, tolerance=tolerance),
               SIMPLIFY=TRUE, USE.NAMES=FALSE)
    }
    attributes(res) <- attributes(x)  # set dim and dimnames
    res
}


### apply a set of scoring rules at once

scores.default <- function(x, mu, size = NULL,
                           which = c("logs", "rps", "dss", "ses"),
                           sign = FALSE, ...)
{
    stopifnot(is.na(mu) | mu >= 0)
    if (!is.null(size)) stopifnot(is.na(size) | size > 0)
    
    ## compute individual scores (these have the same dimensions as x)
    scorelist <- lapply(X = setNames(nm = which), FUN = do.call,
                        args = alist(x = x, mu = mu, size = size),
                        envir = environment())

    ## append sign of x-mu
    if (sign)
        scorelist <- c(scorelist, list("sign" = sign(x-mu)))

    ## gather scores in an array
    simplify2array(scorelist, higher = TRUE)
}


### apply scoring rules to a set of oneStepAhead() forecasts

scores.oneStepAhead <- function (x, which = c("logs","rps","dss","ses"),
                                 units = NULL, sign = FALSE, individual = FALSE,
                                 reverse = FALSE, ...)
{
    y <- x$observed  # observed counts during the prediction window
    mu <- x$pred     # predicted counts (same dim as y)
    ## transform overdispersion to dnbinom() parameterization
    size <- psi2size.oneStepAhead(x) # -> NULL or full dim(y) matrix

    ## select units
    if (!is.null(units)) {
        y <- y[,units,drop=FALSE]
        mu <- mu[,units,drop=FALSE]
        size <- size[,units,drop=FALSE] # works with size = NULL
    }
    nUnits <- ncol(y)
    if (nUnits == 1L)
        individual <- TRUE  # no need to apply rowMeans() below

    result <- scores.default(x = y, mu = mu, size = size,
                             which = which, sign = sign)

    ## reverse order of the time points (historically)
    if (reverse) {
        result <- result[nrow(result):1L,,,drop=FALSE]
    }

    ## average over units if requested
    if (individual) {
        drop(result)
    } else {
        apply(X=result, MARGIN=3L, FUN=rowMeans)
        ## this gives a nrow(y) x (5L+sign) matrix (or a vector in case nrow(y)=1)
    }
}


## calculate scores with respect to fitted values

scores.hhh4 <- function (x, which = c("logs","rps","dss","ses"),
                         subset = x$control$subset, units = seq_len(x$nUnit),
                         sign = FALSE, ...)
{
    ## slow implementation via "fake" oneStepAhead():
    ##fitted <- oneStepAhead(x, tp = subset[1L] - 1L, type = "final",
    ##                       keep.estimates = FALSE, verbose = FALSE)
    ##scores.oneStepAhead(fitted, which = which, units = units, sign = sign,
    ##                    individual = TRUE, reverse = FALSE)

    result <- scores.default(
        x = x$stsObj@observed[subset, units, drop = FALSE],
        mu = x$fitted.values[match(subset, x$control$subset), units, drop = FALSE],
        size = psi2size.hhh4(x, subset, units),
        which = which, sign = sign)
    rownames(result) <- subset
    drop(result)
}
