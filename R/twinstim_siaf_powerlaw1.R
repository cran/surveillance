################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### 1-parameter power-law kernel f(s) = (1 + ||s||)^-d, i.e., sigma = 1
###
### Copyright (C) 2019 Sebastian Meyer
### $Revision: 2430 $
### $Date: 2019-07-02 16:32:45 +0200 (Tue, 02. Jul 2019) $
################################################################################


siaf.powerlaw1 <- function (nTypes = 1, validpars = NULL, sigma = 1)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)
    stopifnot(isScalar(sigma), sigma > 0)
    SIAF <- siaf.powerlaw(nTypes)  # we can reuse some functions from there

    ## for the moment we don't make this type-specific
    if (nTypes != 1) stop("type-specific shapes are not yet implemented")

    ## spatial kernel
    f <- function (s, logd, types = NULL, sigma = 1) {
        d <- exp(logd)
        sLength <- sqrt(.rowSums(s^2, nrow(s), 2L))
        (sLength + sigma)^-d
    }
    ## set desired sigma as default value
    formals(f)$sigma <- sigma
    environment(f) <- baseenv()

    ## numerically integrate f over a polygonal domain
    F <- function (polydomain, f, logd, type = NULL, logsigma = 0, ...) {
        logpars <- c(logsigma, logd)
        siaf_polyCub_iso(polydomain$bdry, "intrfr.powerlaw", logpars,
                         list(...))
    }
    formals(F)$logsigma <- log(sigma)
    environment(F) <- getNamespace("surveillance")

    ## fast integration of f over a circular domain
    Fcircle <- SIAF$Fcircle  # hack original two-parameter version ...
    body(Fcircle)[2:4] <- NULL
    names(formals(Fcircle))[2] <- "logd"
    formals(Fcircle)$sigma <- sigma

    ## derivative of f wrt logpars
    deriv <- function (s, logd, types = NULL, sigma = 1) {
        d <- exp(logd)
        sLength <- sqrt(.rowSums(s^2, nrow(s), 2L))
        tmp <- -d*log(sLength + sigma)
        matrix(tmp * exp(tmp))
    }
    formals(deriv)$sigma <- sigma
    environment(deriv) <- baseenv()

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- function (polydomain, deriv, logd, type = NULL, logsigma = 0, ...) {
        logpars <- c(logsigma, logd)
        siaf_polyCub_iso(polydomain$bdry, "intrfr.powerlaw.dlogd",
                         logpars, list(...))
    }
    formals(Deriv)$logsigma <- log(sigma)
    environment(Deriv) <- getNamespace("surveillance")

    ## Simulation function (via polar coordinates)
    simulate <- SIAF$simulate  # hack original two-parameter version ...
    names(formals(simulate))[2] <- "logd"
    formals(simulate)$logsigma <- log(sigma)
    body(simulate) <- as.call(
        append(as.list(body(simulate)),
               quote(siafpars <- c(logsigma, logd)),
               after = 1)
    )

    ## return the kernel specification
    list(f = f, F = F, Fcircle = Fcircle, deriv = deriv, Deriv = Deriv,
         simulate = simulate, npars = 1L, validpars = validpars)
}
