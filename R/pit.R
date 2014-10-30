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
### $Revision: 1056 $
### $Date: 2014-10-07 12:10:13 +0200 (Tue, 07 Oct 2014) $
################################################################################


## x - observed count data
## pdistr - predictive CDF or a list of such predictive CDF's,
##          one for each data point x. If evaluated at x=-1 it must return 0
## J - number of bins 
## ... - additional arguments for pdistr(), recycled to the length of x.
##       Ignored if pdistr is a list.
## plot - NULL (no plot) or a list of arguments for plot.histogram

pit <- function (x, pdistr, J=10, relative=TRUE, ..., plot = list())
{
    PxPxm1 <- pitPxPxm1(x, pdistr, ...)
    breaks <- (0:J)/J
    Fbar_seq <- vapply(X = breaks, FUN = pit1, FUN.VALUE = 0,
                       Px = PxPxm1[1L,], Pxm1 = PxPxm1[2L,], USE.NAMES = FALSE)
    scale <- if (relative) J else 1
    f_j <- scale * diff.default(Fbar_seq)
    
    res <- structure(list(breaks=breaks, counts=f_j, density=f_j,
                          mids=breaks[-(J+1)] + 1/J/2,
                          xname="PIT", equidist=TRUE),
                     class="histogram")

    if (is.null(plot)) res else pitplot(res, plot)
}

pitPxPxm1 <- function (x, pdistr, ...)
{
    if (is.list(pdistr)) { # list of functions, not necessarily vectorized
        stopifnot(length(pdistr) == length(x))
        vapply(X = seq_along(x), FUN = function (i) {
            stopifnot(isTRUE(
                all.equal.numeric(0, pdistr[[i]](-1), check.attributes = FALSE)
            ))
            c(pdistr[[i]](x[i]), pdistr[[i]](x[i]-1))
        }, FUN.VALUE = c(0,0), USE.NAMES = FALSE) # 2 x length(x)
    } else { # pdistr is (the name of) a function
        pdistr <- match.fun(pdistr)
        if (nargs() == 2L) { # no dots, same pdistr for every data point
                             # and assumed to be vectorized
            stopifnot(isTRUE(all.equal.numeric(0, pdistr(-1))))
            rbind(pdistr(x),
                  pdistr(x-1),
                  deparse.level = 0)
        } else { # ... arguments for pdistr, recycled to the length of x
                 # pdistr is called by mapply, so no need to be vectorized
            stopifnot(isTRUE(all.equal.numeric(
                0, do.call("pdistr", c(list(-1), lapply(list(...), "[", 1L))),
                check.attributes = FALSE)))
            rbind(mapply(pdistr, x, ..., SIMPLIFY = TRUE, USE.NAMES = FALSE),
                  mapply(pdistr, x-1, ..., SIMPLIFY = TRUE, USE.NAMES = FALSE),
                  deparse.level = 0)
        }
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
                        ylab = if(relative) "Relative Frequency" else "Density")
    args <- if (is.list(args)) {
        args[["x"]] <- NULL             # manual x is ignored
        args <- modifyList(defaultArgs, args)
    } else defaultArgs
    do.call("plot", args)
    abline(h=if (relative) 1 else 1/length(pit$mids), lty=2, col="grey")
    invisible(pit)
}
