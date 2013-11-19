################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### _L_agged power-law kernel f(s) = (||s||/sigma)^-d for ||s|| >= sigma, else 1
### Similar to the density of the Pareto distribution (but value 1 for < sigma)
###
### Copyright (C) 2013 Sebastian Meyer
### $Revision: 666 $
### $Date: 2013-11-08 15:45:36 +0100 (Fre, 08 Nov 2013) $
################################################################################


siaf.powerlawL <- function (nTypes = 1, logpars = TRUE,
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
        expression(ifelse(sLength < sigma, 1, (sLength/sigma)^-d))
    ))

    ## numerically integrate f over a polygonal domain
    F <- siaf.fallback.F
    ## alternative (too slow): subtract the inner sigma-region with constant f=1 from polydomain
    ## F <- function (polydomain, f, logpars, type=NULL, nGQ=20L, nCircle2Poly=32)
    ## {
    ##     sigma <- exp(logpars[[1L]])
    ##     intinner <- pi*sigma^2
    ##     poly <- setdiff(owin2gpc(polydomain),
    ##                     discpoly(c(0,0), sigma, npoly=nCircle2Poly, class="gpc.poly"))
    ##     intpoly <- polyCub.SV(poly, f, logpars, type=type, alpha=0, nGQ=nGQ)
    ##     intinner + intpoly
    ## }
    
    ## fast integration of f over a circular domain
    Fcircle <- function (r, logpars, type = NULL) {}
    body(Fcircle) <- as.call(c(as.name("{"),
        tmp,
        expression(
            ## trivial case: radius of integration domain < sigma (=> constant f)
            if (r <= sigma) return(pi * r^2),
            
            ## otherwise, if r > sigma, integration via f^-1
            fofr <- (r/sigma)^-d,
            basevolume <- pi * r^2 * fofr,   # cylinder volume up to height f(r)
            intfinvsq <- sigma^2 * if (d == 2) -d*log(sigma/r) else {
                d/(d-2) * (1 - (sigma/r)^(d-2))
            },
            basevolume + pi * intfinvsq
        )
    ))

    ## derivative of f wrt logpars
    deriv <- function (s, logpars, types = NULL) {}
    body(deriv) <- as.call(c(as.name("{"),
        tmp,
        expression(
            sLength <- sqrt(.rowSums(s^2, nrow(s), 2L)),
            f <- (sLength/sigma)^-d,
            derivlogsigma <- d*f,
            derivlogd <- f*log(f),
            idx0 <- sLength < sigma,
            derivlogsigma[idx0] <- derivlogd[idx0] <- 0,
            cbind(derivlogsigma, derivlogd)
            )
    ))

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- function (polydomain, deriv, logpars, type = NULL, nGQ = 20L) {}
    body(Deriv) <- as.call(c(as.name("{"),
        ## Determine a = argmax(abs(deriv(c(x,0))))
        c(tmp, expression(
            a.logsigma <- sigma,
            a.logd <- sigma * exp(1/d),
            xrange <- polydomain$xrange,           # polydomain is a "owin"
            maxxdist <- max(abs(xrange)),
            a <- pmin(c(a.logsigma, a.logd), maxxdist),
            if (sum(xrange) < 0) a <- -a # is more of the domain left of 0?
            )),
        expression(
            deriv1 <- function (s, paridx)
                deriv(s, logpars, type)[,paridx,drop=TRUE],
            intderiv1 <- function (paridx)
                polyCub.SV(polydomain, deriv1, paridx=paridx,
                           nGQ = nGQ, alpha = a[paridx]),
            res.logsigma <- intderiv1(1L),
            res.logd <- intderiv1(2L),
            res <- c(res.logsigma, res.logd),
            res
            )
    ))

    ## "effective" integration range (based on quantile of the Pareto distri)
    effRange <- if (isScalar(effRangeProb)) {
        stop("'effRange' is currently not supported for power law's")
        effRange <- function (logpars) {}
        body(effRange) <- as.call(c(as.name("{"),
            tmp,
            expression(
                alpha <- d-1,  # only works for alpha > 0, i.e. d > 1 !
                sigma/(1-effRangeProb)^(1/alpha)
                )
        ))
        effRange
    } else NULL

    ## simulate from the lagged power-law kernel (within a maximum distance 'ub')
    simulate <- function (n, logpars, type, ub)
    {
        sigma <- exp(logpars[[1L]])
        d <- exp(logpars[[2L]])
        ## Sampling via polar coordinates and inversion method

        ## random angle
        theta <- stats::runif(n, 0, 2*pi)

        ## sampling radius r
        ## trivial case u < sigma: p(r) \propto r on [0;u]
        if (ub < sigma) {
            r <- ub * sqrt(stats::runif(n)) # inversion sampling
            ## now rotate each point by a random angle to cover all directions
            return(r * cbind(cos(theta), sin(theta)))
        }

        ## case u >= sigma: p(r) \propto r if r<sigma, r*(r/sigma)^-d otherwise
        ## sample hierarchically from mixture
        ## calculate probability for r < sigma (uniform short-range component)
        mass1 <- sigma^2/2              # = int_0^sigma x dx
        ## mass2 = int_sigma^u x * (x/sigma)^-d dx; corresponding primitive:
        ## prim <- function (x) {
        ##     sigma^d * if (d == 2) log(x) else x^(2-d) / (2-d)
        ## }
        mass2 <- sigma^d *
            if (d == 2) log(ub/sigma) else (ub^(2-d)-sigma^(2-d))/(2-d)
        ## probability for r < sigma is mass1/(mass1+mass2) => sample component
        unir <- stats::runif(n) <= mass1 / (mass1 + mass2)

        ## samples from the uniform short-range component:
        n1 <- sum(unir)
        r1 <- sigma * sqrt(stats::runif(n1)) # similar to the case u < sigma

        ## samples from power-law component: p2(r) \propto r^(-d+1) on [sigma;u]
        ## For d>2 only, we could use VGAM::rpareto(n,sigma,d-2), d=1 is trivial
        n2 <- n - n1
        r2 <- if (d==1) stats::runif(n2, sigma, ub) else { # inversion sampling
            P2inv <- if (d == 2) { function (z) ub^z * sigma^(1-z) } else {
                function (z) (z*ub^(2-d) + (1-z)*sigma^(2-d))^(1/(2-d))
            }
            P2inv(stats::runif(n2))
        }

        ## put samples from both components together
        r <- c(r1, r2)
        
        ## now rotate each point by a random angle to cover all directions
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
