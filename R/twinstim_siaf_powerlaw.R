################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Power-law kernel f(s) = (||s||+sigma)^-d
### This is the pure kernel of the Lomax density (the density requires d>1, but
### for the siaf specification we only want d to be positive)
###
### Copyright (C) 2013-2014,2017 Sebastian Meyer
### $Revision: 1988 $
### $Date: 2017-10-06 11:04:19 +0200 (Fri, 06. Oct 2017) $
################################################################################


siaf.powerlaw <- function (nTypes = 1, validpars = NULL, engine = "C")
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)
    engine <- match.arg(engine, c("C", "R"))

    ## for the moment we don't make this type-specific
    if (nTypes != 1) stop("type-specific shapes are not yet implemented")

    ## helper expression, note: logpars=c(logscale=logsigma, logd=logd)
    tmp <- expression(
        logsigma <- logpars[[1L]],  # used "[[" to drop names
        logd <- logpars[[2L]],
        sigma <- exp(logsigma),
        d <- exp(logd)
        )

    ## spatial kernel
    f <- function (s, logpars, types = NULL) {}
    body(f) <- as.call(c(as.name("{"),
        tmp,
        expression(sLength <- sqrt(.rowSums(s^2, nrow(s), 2L))),
        expression((sLength+sigma)^-d)
    ))
    environment(f) <- baseenv()

    ## numerically integrate f over a polygonal domain
    F <- siaf_F_polyCub_iso(intrfr_name = "intrfr.powerlaw", engine = engine)

    ## fast integration of f over a circular domain
    Fcircle <- function (r, logpars, type = NULL) {}
    body(Fcircle) <- as.call(c(as.name("{"),
        tmp,
        expression(
            fofr <- (r+sigma)^-d,
            fof0 <- sigma^-d,
            ## calculate cylinder volume up to height f(r)
            basevolume <- if (is.infinite(r)) 0 else pi * r^2 * fofr,
            ## r=Inf is used in R0(,trimmed=F), Fcircle(Inf) is finite if d>2
            Ifinvsq <- function (z) {
                if (d == 1) {
                    -1/z - 2*sigma*log(z) + sigma^2*z
                } else if (d == 2) {
                    log(z) - 4*sigma*sqrt(z) + sigma^2*z
                } else {
                    z^(1-2/d) * d / (d-2) - z^(1-1/d) * 2*sigma*d/(d-1) + sigma^2*z
                }
            },
            intfinvsq <- Ifinvsq(fof0) - Ifinvsq(fofr),
            basevolume + pi * intfinvsq
            )
    ))
    environment(Fcircle) <- baseenv()

    ## derivative of f wrt logpars
    deriv <- function (s, logpars, types = NULL) {}
    body(deriv) <- as.call(c(as.name("{"),
        tmp,
        expression(
            sLength <- sqrt(.rowSums(s^2, nrow(s), 2L)),
            rsigma <- sLength + sigma,
            rsigmad <- rsigma^d,
            derivlogsigma <- -d*sigma / rsigmad / rsigma,
            derivlogd <- -d*log(rsigma) / rsigmad,
            cbind(derivlogsigma, derivlogd)
            )
    ))
    environment(deriv) <- baseenv()

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- siaf_Deriv_polyCub_iso(
        intrfr_names = c("intrfr.powerlaw.dlogsigma", "intrfr.powerlaw.dlogd"),
        engine = engine)

    ## Simulation function (via polar coordinates)
    simulate <- siaf.simulatePC(intrfr.powerlaw)
    ## if (!is.finite(ub)) normconst <- {
    ##     ## for sampling on [0;Inf] the density is only proper if d > 2
    ##     if (d <= 2) stop("improper density for d<=2, 'ub' must be finite")
    ##     1/(sigma^(d-2) * (d-2)*(d-1)) # = intrfr.powerlaw(Inf)
    ## }
    environment(simulate) <- getNamespace("surveillance")

    ## return the kernel specification
    list(f=f, F=F, Fcircle=Fcircle, deriv=deriv, Deriv=Deriv,
         simulate=simulate, npars=2L, validpars=validpars)
}


## integrate x*f(x) from 0 to R (vectorized)
intrfr.powerlaw <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    if (d == 1) {
        R - sigma * log(R/sigma + 1)
    } else if (d == 2) {
        log(R/sigma + 1) - R/(R+sigma)
    } else {
        (R*(R+sigma)^(1-d) - ((R+sigma)^(2-d) - sigma^(2-d))/(2-d)) / (1-d)
    }
}
## local({ # validation via numerical integration -> tests/testthat/test-siafs.R
##     p <- function (r, sigma, d) r * (r+sigma)^-d
##     Pnum <- function (r, sigma, d) sapply(r, function (.r) {
##         integrate(p, 0, .r, sigma=sigma, d=d)$value
##     })
##     r <- c(1,2,5,10,20,50,100)
##     dev.null <- sapply(c(1,2,1.6), function(d) stopifnot(isTRUE(
##         all.equal(intrfr.powerlaw(r, log(c(3, d))), Pnum(r, 3, d)))))
## })

## integrate x * (df(x)/dlogsigma) from 0 to R (vectorized)
intrfr.powerlaw.dlogsigma <- function (R, logpars, types = NULL)
{
    pars <- exp(logpars)
    -prod(pars) * intrfr.powerlaw(R, log(pars+c(0,1)), types)
}

## integrate x * (df(x)/dlogd) from 0 to R (vectorized)
## (thanks to Maple 17) -> validated in tests/testthat/test-siafs.R
intrfr.powerlaw.dlogd <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    if (d == 1) {
        sigma * logpars[[1L]] * (1-logpars[[1L]]/2) - log(R+sigma) * (R+sigma) +
            sigma/2 * log(R+sigma)^2 + R
    } else if (d == 2) {
        (-log(R+sigma) * ((R+sigma)*log(R+sigma) + 2*sigma) +
         (R+sigma)*logpars[[1L]]*(logpars[[1L]]+2) + 2*R) / (R+sigma)
    } else {
        (sigma^(2-d) * (logpars[[1L]]*(-d^2 + 3*d - 2) - 2*d + 3) +
         (R+sigma)^(1-d) * (log(R+sigma)*(d-1)*(d-2) * (R*(d-1) + sigma) +
                            R*(d^2+1) + 2*d*(sigma-R) - 3*sigma)
         ) * d / (d-1)^2 / (d-2)^2
    }
}
