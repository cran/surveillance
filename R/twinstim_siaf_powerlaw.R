################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Power-law kernel f(s) = (||s||+sigma)^-d
### This is the pure kernel of the Lomax density (the density requires d>1, but
### for the siaf specification we only want d to be positive)
###
### Copyright (C) 2013 Sebastian Meyer
### $Revision: 631 $
### $Date: 2013-08-28 17:52:52 +0200 (Mit, 28 Aug 2013) $
################################################################################


siaf.powerlaw <- function (nTypes = 1, logpars = TRUE,
                           effRangeProb = NULL, validpars = NULL)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)

    ## for the moment we don't make this type-specific
    if (nTypes != 1) stop("type-specific shapes are not yet implemented")
    if (!logpars) stop("only the 'logpars' parametrization is implemented")

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

    ## numerically integrate f over a polygonal domain
    F <- siaf.fallback.F
    
    ## fast integration of f over a circular domain
    Fcircle <- function (r, logpars, type = NULL) {}
    body(Fcircle) <- as.call(c(as.name("{"),
        tmp,
        expression(
            fofr <- (r+sigma)^-d,
            fof0 <- sigma^-d,
            basevolume <- pi * r^2 * fofr, # cylinder volume up to height f(r)
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

    ## derivative of f wrt logpars
    deriv <- function (s, logpars, types = NULL) {}
    body(deriv) <- as.call(c(as.name("{"),
        tmp,
        expression(
            sLength <- sqrt(.rowSums(s^2, nrow(s), 2L)),
            rsigmad <- (sLength+sigma)^d,
            derivlogsigma <- -d*sigma / rsigmad / (sLength+sigma),
            derivlogd <- -log(rsigmad) / rsigmad,
            cbind(derivlogsigma, derivlogd)
            )
    ))

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- function (polydomain, deriv, logpars, type = NULL, nGQ = 20L) {}
    body(Deriv) <- as.call(c(as.name("{"),
        ## Determine a = argmax(abs(deriv(c(x,0))))
        c(tmp, expression(
            xrange <- polydomain$xrange,           # polydomain is a "owin"
            maxxdist <- max(abs(xrange)),
            a.logsigma <- 0,
            xmin.derivlogd <- exp(1/d) - sigma,
            a.logd <- if (xmin.derivlogd <= 0) 0 else {
                derivlogd.0 <- -d*log(sigma) / sigma^d
                derivlogd.xmin <- -exp(-1)
                if (abs(derivlogd.0) >= abs(derivlogd.xmin))
                    0 else min(maxxdist, xmin.derivlogd)
            },
            a <- c(a.logsigma, a.logd),
            if (sum(xrange) < 0) a <- -a # is more of the domain left of 0?
            )),
        expression(
            deriv1 <- function (s, paridx)
                deriv(s, logpars, type)[,paridx,drop=TRUE],
            intderiv1 <- function (paridx)
                polyCub::polyCub.SV(polydomain, deriv1, paridx=paridx,
                           nGQ = nGQ, alpha = a[paridx]),
            res.logsigma <- intderiv1(1L),
            res.logd <- intderiv1(2L),
            res <- c(res.logsigma, res.logd),
            res
            )
    ))

    ## "effective" integration range (based on some high quantile)
    effRange <- if (isScalar(effRangeProb)) {
        stop("'effRange' is currently not supported for power law's")
        effRange <- function (logpars) {}
        body(effRange) <- as.call(c(as.name("{"),
            substitute(qlomax(effRangeProb, exp(logpars[[1]]), exp(logpars[[2]])-1),
                       list(effRangeProb=effRangeProb)) # only works for d > 1!
        ))
        effRange
    } else NULL

    ## simulate from the power-law kernel (within a maximum distance 'ub')
    simulate <- function (n, logpars, type, ub)
    {
        sigma <- exp(logpars[[1L]])
        d <- exp(logpars[[2L]])

        ## Sampling from f_{2D}(s) \propto f(||s||), truncated to ||s|| <= ub,
        ## works via polar coordinates and the inverse transformation method:
        ## p_{2D}(r,theta) = r * f_{2D}(x,y) \propto r*f(r) = r*(r+\sigma)^{-d}
        ## => sample angle theta from U(0,2*pi) and r according to:
        ## p <- function (r) r * (r+sigma)^-d
        
        ## Primitive of p: P(r) = int_0^r p(z) dz
        P <- function (r) {
            if (d == 1) {
                r - sigma * log(r/sigma + 1)
            } else if (d == 2) {
                log(r/sigma + 1) - r/(r+sigma)
            } else {
                (r*(r+sigma)^(1-d) - (r+sigma)^(2-d) / (2-d) +
                 sigma^(2-d)/(2-d)) / (1-d)
            }
        }
        ## for validation of analytic P: numerical integration
        ## Pnum <- function (r, sigma, d) sapply(r, function (.r) {
        ##     integrate(p, 0, .r, sigma=sigma, d=d)$value
        ## })        

        ## Normalizing constant of the density on [0;ub]
        normconst <- if (is.finite(ub)) P(ub) else {
            ## for sampling on [0;Inf] the density is only proper if d > 2
            if (d <= 2) stop("improper density for d<=2, 'ub' must be finite")
            (sigma^(d-2) * (d-2) * (d-1))^-1 # = P(Inf)
        }

        ## => cumulative distribution function
        CDF <- function (q) P(q) / normconst

        ## For inversion sampling, we need the quantile function CDF^-1
        ## However, this is not available in closed form, so we use uniroot,
        ## which requires a finite upper bound!
        ## Note: in simEpidataCS, simulation is always bounded to eps.s and to
        ## the largest extend of W, thus, 'ub' is finite
        stopifnot(is.finite(ub))
        QF <- function (p) uniroot(function(q) CDF(q)-p, lower=0, upper=ub)$root

        ## Now sample r as QF(U), where U ~ U(0,1)
        r <- sapply(runif(n), QF)
        ## Check simulation via kernel estimate:
        ## plot(density(r, from=0, to=ub)); curve(p(x) / normconst, add=TRUE, col=2)
        
        ## now rotate each point by a random angle to cover all directions
        theta <- stats::runif(n, 0, 2*pi)
        r * cbind(cos(theta), sin(theta))
    }

    ## set function environments to the global environment
    environment(f) <- environment(F) <- environment(Fcircle) <-
        environment(deriv) <- environment(Deriv) <-
            environment(simulate) <- .GlobalEnv
    if (is.function(effRange)) environment(effRange) <- .GlobalEnv

    ## return the kernel specification
    list(f=f, F=F, Fcircle=Fcircle, effRange=effRange, deriv=deriv, Deriv=Deriv,
         simulate=simulate, npars=2L, validpars=validpars)
}
