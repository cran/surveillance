################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Non-randomized version of the PIT histogram as discussed in:
### Predictive model assessment for count data
### Czado, C., Gneiting, T. & Held, L. (2009)
### Biometrics 65:1254-1261
###
### Copyright (C) 2010-2012 Michaela Paul, 2013-2014 Sebastian Meyer
### $Revision: 930 $
### $Date: 2014-05-22 17:16:59 +0200 (Thu, 22 May 2014) $
################################################################################


## x - observed data
## pdistr - predictive CDF, i.e. a vectorized function (x, ...)
##          or a list of such predictive CDF's, one for each data point x
##          If evaluated at x=-1 it should return 0
## J - number of bins 
## ... - arguments for pdistr. Ignored if pdistr is a list.
## plot - NULL (no plot) or a list of arguments for plot.histogram

pit <- function (x, pdistr, J=10, relative=TRUE, ..., plot = list())
{
    PxPxm1 <- pitPxPxm1(x, pdistr, ...)
    breaks <- (0:J)/J
    Fbar_seq <- sapply(breaks, pit1, Px=PxPxm1[1,], Pxm1=PxPxm1[2,])
    scale <- if (relative) J else 1
    f_j <- scale * diff(Fbar_seq)
    
    res <- structure(list(breaks=breaks, counts=f_j, density=f_j,
                          mids=breaks[-(J+1)] + 1/J/2,
                          xname="PIT", equidist=TRUE),
                     class="histogram")

    if (is.null(plot)) res else pitplot(res, plot)
}

pitPxPxm1 <- function (x, pdistr, ...)
{
    if (is.list(pdistr)) {
        stopifnot(length(pdistr) == length(x))
        sapply(seq_along(x), function (i) {
            stopifnot(isTRUE(all.equal(0, pdistr[[i]](-1))))
            res <- pdistr[[i]](c(x[i], x[i]-1))
            if (length(res) == 2 && is.vector(res, mode="numeric")) res else
                stop("pdistr[[", i, "]](c(", x[i], ", ", x[i]-1,
                     ")) is no numeric vector of length 2")
        }) # 2 x length(x)
    } else { # same pdistr for every data point
        stopifnot(pdistr(-1, ...) == 0)
        rbind(pdistr(x, ...), pdistr(x-1, ...))
    }
}

## calculate \bar{F}(u) for scalar u
pit1 <- function (u, Px, Pxm1)
{
    if (u <= 0) return(0) else if (u >= 1) return(1)
    F_u <- (u-Pxm1) / (Px-Pxm1)
    ## If Px=Pxm1, this means that predict. prob. of observed x is exactly zero.
    ## We get NaN for F_u. Our predictive model is bad if that happens.
    ## We could assign either 0 or 1 to express that and issue a warning.
    if (any(is.nan(F_u))) {
        warning("predictive distribution has 0 probability for observed 'x'")
        F_u[is.nan(F_u)] <- 0
    }
    F_u[F_u < 0] <- 0
    F_u[F_u > 1] <- 1
    mean(F_u)
}

## plot the PIT histogram
pitplot <- function (pit, args)
{
    relative <- !isTRUE(all.equal(1, sum(pit$density)))
    defaultArgs <- list(x = pit, main = "",
                        ylab = if(relative) "Relative frequency" else "Density")
    args <- if (is.list(args)) {
        args[["x"]] <- NULL             # manual x is ignored
        args <- modifyList(defaultArgs, args)
    } else defaultArgs
    do.call("plot", args)
    abline(h=if (relative) 1 else 1/length(pit$mids), lty=2, col="grey")
    invisible(pit)
}
