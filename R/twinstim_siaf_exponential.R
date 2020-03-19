################################################################################
### Exponential kernel f(s) = exp(-||s||/sigma)
###
### Copyright (C) 2020 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################


siaf.exponential <- function (nTypes = 1, validpars = NULL, engine = "C")
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)
    engine <- match.arg(engine, c("C", "R"))

    ## for the moment we don't make this type-specific
    if (nTypes != 1) stop("type-specific shapes are not yet implemented")

    ## spatial kernel
    f <- function (s, logsigma, types = NULL) {
        sigma <- exp(logsigma)
        sLength <- sqrt(.rowSums(s^2, nrow(s), 2L))
        exp(-sLength/sigma)
    }
    environment(f) <- baseenv()

    ## numerically integrate f over a polygonal domain
    F <- siaf_F_polyCub_iso(intrfr_name = "intrfr.exponential", engine = engine)

    ## fast integration of f over a circular domain
    Fcircle <- function (r, logsigma, type = NULL) {
        sigma <- exp(logsigma)
        fofr <- exp(-r/sigma)
        ## f(r) approaches 0, r=Inf is used in R0(,trimmed=F)
        if (fofr == 0)
            return(2*pi*sigma^2)
        ## calculate cylinder volume up to height f(r)
        basevolume <- pi * r^2 * fofr
        ## integration via f^-1
        Ifinvsq <- function (z) sigma^2 * z * ((log(z)-1)^2 + 1)
        ##intfinvsq <- Ifinvsq(fof0) - Ifinvsq(fofr)  # fof0 = 1
        intfinvsq <- 2*sigma^2 - Ifinvsq(fofr)
        basevolume + pi * intfinvsq
    }
    environment(Fcircle) <- baseenv()

    ## derivative of f wrt logsigma
    deriv <- function (s, logsigma, types = NULL) {
        sigma <- exp(logsigma)
        sLength <- sqrt(.rowSums(s^2, nrow(s), 2L))
        z <- sLength/sigma
        matrix(z * exp(-z))
    }
    environment(deriv) <- baseenv()

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- siaf_Deriv_polyCub_iso(intrfr_names = "intrfr.exponential.dlogsigma",
                                    engine = engine)

    ## simulation from the kernel (via polar coordinates)
    simulate <- siaf.simulatePC(intrfr.exponential)
    environment(simulate) <- getNamespace("surveillance")

    ## return the kernel specification
    list(f = f, F = F, Fcircle = Fcircle, deriv = deriv, Deriv = Deriv,
         simulate = simulate, npars = 1L, validpars = validpars)
}


## integrate x*f(x) from 0 to R (vectorized)
intrfr.exponential <- function (R, logsigma, types = NULL)
{
    sigma <- exp(logsigma)
    sigma * (sigma - (R+sigma)*exp(-R/sigma))
}

## integrate x * (df(x)/dlogsigma) from 0 to R (vectorized)
## Note: df(x)/dlogsigma = x * exp(-(x/sigma)-logsigma)
intrfr.exponential.dlogsigma <- function (R, logsigma, types = NULL)
{
    sigma <- exp(logsigma)
    2*sigma^2 - ((R+sigma)^2 + sigma^2)*exp(-R/sigma)
}
