################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Implementation of an isotropic power-law kernel:
### Density of the Lomax distribution, i.e. a shifted Pareto distribution with
### domain [0;Inf), scaled to 1 at the origin.
### THIS IS OBSOLETE since siaf.powerlaw provides the same kernel but with a
### better parametrization (using d := alpha+1 and omitting prop. constants)
###
### Copyright (C) 2012-2013 Sebastian Meyer
### $Revision: 527 $
### $Date: 2013-03-08 13:25:48 +0100 (Fr, 08. Mrz 2013) $
################################################################################


## density=FALSE returns standardized Lomax kernel, i.e. f(x) = f_Lomax(x)/f(0),
## such that the kernel function starts at 1. f_Lomax(0) = alpha / sigma

siaf.lomax <- function (nTypes = 1, logpars = TRUE, density = FALSE,
                        effRangeProb = NULL, validpars = NULL)
{
    .Deprecated("siaf.powerlaw")
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)

    ## for the moment we don't make this type-specific
    if (nTypes != 1) stop("type-specific shapes are not yet implemented")
    if (!logpars) stop("only the 'logpars' parametrization is implemented")

    ## helper expression, note: logpars=c(logscale=logsigma, logshape=logalpha)
    tmp <- expression(
        logsigma <- logpars[[1L]],  # used "[[" to drop names
        logalpha <- logpars[[2L]],
        sigma <- exp(logsigma),
        alpha <- exp(logalpha)
        )

    ## spatial kernel
    f <- function (s, logpars, types = NULL) {}
    body(f) <- as.call(c(as.name("{"),
        tmp,
        expression(sLength <- sqrt(rowSums(s^2))),
        expression(logfvals <- (alpha+1) * (logsigma - log(sLength+sigma))),
        ##logfvals <- logalpha + alpha*logsigma - (alpha+1)*log(sLength+sigma)
        if (density) {
            expression(exp(logfvals + logalpha - logsigma))
        } else {
            expression(exp(logfvals))
        }
    ))

    ## numerically integrate f over a polygonal domain
    F <- siaf.fallback.F
    
    ## fast integration of f over a circular domain
    Fcircle <- function (r, logpars, type = NULL) {}
    body(Fcircle) <- as.call(c(as.name("{"),
        tmp,
        ## integrate standardized f(x) = f_Lomax(x) / f_Lomax(0)
        expression(
            logfofr <- (alpha+1) * (logsigma - log(r+sigma)),
            fofr <- exp(logfofr),
            basevolume <- pi * r^2 * fofr, # cylinder volume up to height f(r)
            intfinvsq <- if (alpha == 1) {
                sigma^2 * (-3 - logfofr + 4*sqrt(fofr) - fofr)
            } else {
                (2*sigma^2 *(1-fofr) - r*(alpha+1)*(alpha*r+2*sigma) *fofr) /
                    alpha / (alpha-1)
            },
            int <- basevolume + pi * intfinvsq
            ),
         ## if density=TRUE, multiply result by f_Lomax(0)
         if (density) expression(int * alpha / sigma) else expression(int)
    ))

    ## derivative of f wrt logpars
    deriv <- function (s, logpars, types = NULL) {}
    body(deriv) <- as.call(c(as.name("{"),
        tmp,
        expression(sLength <- sqrt(rowSums(s^2))),
        expression(logsigma.xsigma <- logsigma - log(sLength+sigma)),
        if (density) {
            expression(
                logfac <- logalpha - logsigma + (alpha+1) * logsigma.xsigma,
                derivlogsigma.part <- (alpha*sLength-sigma) / (sLength+sigma),
                derivlogalpha.part <- 1 + alpha * logsigma.xsigma
                )
        } else {
            expression(
                logfac <- (alpha+1) * logsigma.xsigma,
                derivlogsigma.part <- (alpha+1) * sLength / (sLength+sigma),
                derivlogalpha.part <- alpha * logsigma.xsigma
                )
        },
        expression(exp(logfac) * cbind(derivlogsigma.part, derivlogalpha.part))
    ))

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- function (polydomain, deriv, logpars, type = NULL, nGQ = 20L) {}
    body(Deriv) <- as.call(c(as.name("{"),
        ## Determine a = argmax(abs(deriv(c(x,0))))
        if (density) {
            expression(a <- c(0,0)) # max abs values are at 0
        } else c(tmp, expression(
            xrange <- polydomain$xrange,           # polydomain is a "owin"
            maxxdist <- max(abs(xrange)),
            a.logsigma <- min(maxxdist, sigma / (alpha+1)),
            a.logalpha <- min(maxxdist, sigma * (exp(1/(alpha+1)) - 1)),
            a <- c(a.logsigma, a.logalpha),
            if (sum(xrange) < 0) a <- -a # is more of the domain left of 0?
            )),
        expression(
            deriv1 <- function (s, paridx)
                deriv(s, logpars, type)[,paridx,drop=TRUE],
            intderiv1 <- function (paridx)
                polyCub.SV(polydomain, deriv1, paridx=paridx,
                           nGQ = nGQ, alpha = a[paridx]),
            res.logsigma <- intderiv1(1L),
            res.logalpha <- intderiv1(2L),
            res <- c(res.logsigma, res.logalpha),
            res
            )
    ))

    ## "effective" integration range (based on some high quantile)
    effRange <- if (isScalar(effRangeProb)) {
        effRange <- function (logpars) {}
        body(effRange) <- as.call(c(as.name("{"),
            substitute(qlomax(effRangeProb, exp(logpars[[1]]), exp(logpars[[2]])),
                       list(effRangeProb=effRangeProb))
        ))
        effRange
    } else NULL

    ## simulate from the Lomax kernel (within a maximum distance 'ub')
    simulate <- function (n, logpars, type, ub)
    {
        ## Sampling from f(s) = dlomax(||s||), truncated to ||s|| <= ub is via
        ## polar coordinates and applies the inverse transformation method
        sigma <- exp(logpars[[1L]])
        alpha <- exp(logpars[[2L]])
        E <- sigma / (alpha-1)        # only valid if alpha > 1
        ##r <- VGAM::rlomax(n, scale=sigma, shape3.q=alpha)
        ## NO, r must be sampled from a density proportional to r*dlomax(r)!!!
        ##curve(x * alpha/sigma * (1+x/sigma)^-(alpha+1), 0, 200)

        ## We need the primitive of the density, i.e. int_0^r x*dlomax(x) dx
        cumd <- if (alpha == 1) {
            function (r) sigma * (log(1+r/sigma) - r/(r+sigma))
        } else {
            function (r)
                ifelse(is.infinite(r) & alpha > 1, E,
                       E * (1 - (1+r/sigma)^-alpha * (r*alpha/sigma+1)))
        }

        ## cumulative distribution function: int_0^q x*dlomax(x) / c dx
        ## the normalization constant c is only finite if alpha > 1 _or_ if 
        ## the density is truncated (in simEpidataCS, simulation is always
        ## bounded to eps.s and to the largest extend of W), otherwise the
        ## distribution is improper
        CDF <- if (is.finite(ub)) { # in this case any alpha > 0 is fine
            normconst <- cumd(ub)
            function (q) cumd(q) / normconst
        } else { # for r in [0;Inf] the density is only proper if alpha > 1
            stopifnot(alpha > 1)
            function (q) cumd(q) / E
        }
        
        ## However, there is no closed form expression for the quantile function
        ## of that distribution (inverse CDF), so we have to use uniroot
        ## However, uniroot needs a finite upper bound
        stopifnot(is.finite(ub))
        QF <- function(p) uniroot(function(q) CDF(q)-p, lower=0, upper=ub)$root

        ## now sample r as QF(U), where U ~ U(0,1)
        r <- sapply(runif(n), QF)
        ## Check simulation via kernel estimate:
        ## plot(density(r, from=0, to=ub))
        ## curve(x*VGAM::dlomax(x,sigma,alpha) / normconst, add=TRUE, col=2)
        
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
