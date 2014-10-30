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
### Copyright (C) 2010-2012 Michaela Paul, 2014 Sebastian Meyer
### $Revision: 1043 $
### $Date: 2014-10-03 09:46:12 +0200 (Fri, 03 Oct 2014) $
################################################################################


## logarithmic score
## logs(P,x) = -log(P(X=x))

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

dss <- function (x, mu, size=NULL) {
    sigma2 <- if (is.null(size)) mu else mu * (1 + mu/size)
    ((x-mu)^2)/sigma2 + log(sigma2)
}


## ranked probability score
## rps(P,x) = sum_0^Kmax {P(X<=k) - 1(x<=k)}^2

rps.one <- function (x, mu, size=NULL, k=40, tolerance=sqrt(.Machine$double.eps)) {
    ## return NA for non-convergent fits (where mu=NA)
    if (is.na(mu)) return(NA_real_)

    ## determine variance of distribution
    sigma2 <- if (is.null(size)) mu else mu * (1 + mu/size)
    
    ## determine the maximum number of summands as Kmax=mean+k*sd
    kmax <- ceiling(mu + k*sqrt(sigma2))
    
    ## compute 1(x<=k)
    ind <- 1*(x < seq_len(kmax+1))
	
    ## compute P(X<=k)
    Px <- if (is.null(size)) {
        ppois(0:kmax, lambda=mu)
    } else {
        pnbinom(0:kmax, mu=mu, size=size)
    }
	
    ## check precision
    if ((1-tail(Px,1))^2 > tolerance)
        warning("precision of finite sum not smaller than tolerance=", tolerance)
	
    ## compute rps
    sum((Px-ind)^2)
}

rps <- function (x, mu, size=NULL, k=40, tolerance=sqrt(.Machine$double.eps)) {
    res <- if (is.null(size)) {
        mapply(rps.one, x=x, mu=mu,
               MoreArgs=list(k=k, tolerance=tolerance), SIMPLIFY=TRUE, USE.NAMES=FALSE)
    } else {
        mapply(rps.one, x=x, mu=mu, size=size,
               MoreArgs=list(k=k, tolerance=tolerance), SIMPLIFY=TRUE, USE.NAMES=FALSE)
    }
    attributes(res) <- attributes(x)  # set dim and dimnames
    res
}


## returns scores in reversed (!) order, i.e. for time points n, n-1, n-2, ...

scores <- function (object, which = c("logs","rps","dss","ses"), units = NULL,
                    sign = FALSE, individual = FALSE)
{
    mu <- object$pred     # predicted counts
    x <- object$observed  # observed counts
    size <- object$psi    # estimated -log(overdispersion), 1 or ncol(x) columns
    ntps <- nrow(x)       # the number of predicted time points
    
    if (!is.null(size)) { # => NegBin
        size <- exp(size) # transform to parameterization suitable for dnbinom()
        if (ncol(size) != ncol(x)) { # => ncol(size)=1, unit-independent psi
            ## replicate to obtain a ntps x ncol(x) matrix
            size <- matrix(size, nrow=ntps, ncol=ncol(x), byrow=FALSE)
        }
        colnames(size) <- colnames(x)  # such that we can select by unit name
    }
    ## At this point, mu, x and size all are ntps x ncol(x) matrices

    ## select units
    if (!is.null(units)) {
        x <- x[,units,drop=FALSE]
        mu <- mu[,units,drop=FALSE]
        size <- size[,units,drop=FALSE]
    }
    nUnits <- ncol(x)
    if (nUnits == 1L)
        individual <- TRUE  # no need to apply rowMeans() below

    ## compute sign of x-mu
    signXmMu <- if(sign) sign(x-mu) else NULL
    
    ## compute individual scores (these are ntps x nUnits matrices)
    scorelist <- lapply(which, do.call, args = alist(x=x, mu=mu, size=size),
                        envir = environment())
    
    ## gather individual scores in an array
    result <- array(c(unlist(scorelist, recursive=FALSE, use.names=FALSE),
                      signXmMu),
                    dim = c(ntps, nUnits, length(which) + sign),
                    dimnames = c(dimnames(x),
                                 list(c(which, if (sign) "sign"))))

    ## reverse order of the time points (historically)
    result <- result[ntps:1L,,,drop=FALSE]

    ## average over units if requested
    if (individual) {
        drop(result)
    } else {
        apply(X=result, MARGIN=3L, FUN=rowMeans)
        ## this gives a ntps x (5L+sign) matrix (or a vector in case ntps=1)
    }
}
