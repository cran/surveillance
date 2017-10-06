################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### _L_agged power-law kernel f(s) = (||s||/sigma)^-d for ||s|| >= sigma, else 1
### Similar to the density of the Pareto distribution (but value 1 for < sigma)
###
### Copyright (C) 2013-2014,2017 Sebastian Meyer
### $Revision: 1988 $
### $Date: 2017-10-06 11:04:19 +0200 (Fri, 06. Oct 2017) $
################################################################################


siaf.powerlawL <- function (nTypes = 1, validpars = NULL, engine = "C")
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
    body(f) <- as.call(c(as.name("{"), tmp, expression(
        sLength <- sqrt(.rowSums(s^2, L <- length(s)/2, 2L)),
        fvals <- rep.int(1, L),
        inPLrange <- which(sLength > sigma),
        fvals[inPLrange] <- (sLength[inPLrange]/sigma)^-d,
        fvals
        )))
    environment(f) <- baseenv()

    ## numerically integrate f over a polygonal domain
    F <- siaf_F_polyCub_iso(intrfr_name = "intrfr.powerlawL", engine = engine)

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
    environment(Fcircle) <- baseenv()

    ## derivative of f wrt logpars
    ## CAVE: the derivative of f wrt logsigma is mathematically NaN at x=sigma
    ## this non-differentiability at the treshhold causes false convergence
    ## warnings by nlminb but is otherwise not relevant (could use slow and
    ## robust Nelder-Mead instead)
    deriv <- function (s, logpars, types = NULL) {}
    body(deriv) <- as.call(c(as.name("{"), tmp, expression(
        sLength <- sqrt(.rowSums(s^2, L <- length(s)/2, 2L)),
        derivlogsigma <- derivlogd <- numeric(L),
        inPLrange <- which(sLength > sigma),
        fPL <- (sLength[inPLrange]/sigma)^-d,
        derivlogsigma[inPLrange] <- d * fPL,
        derivlogd[inPLrange] <- fPL * log(fPL),
        cbind(derivlogsigma, derivlogd)
        )))
    environment(deriv) <- baseenv()

    ## Numerical integration of 'deriv' over a polygonal domain
    Deriv <- siaf_Deriv_polyCub_iso(
        intrfr_names = c("intrfr.powerlawL.dlogsigma", "intrfr.powerlawL.dlogd"),
        engine = engine)

    ## simulate from the lagged power law (within a maximum distance 'ub')
    ##simulate <- siaf.simulatePC(intrfr.powerlawL) # <- generic simulator
    ##environment(simulate) <- getNamespace("surveillance")
    ## faster implementation taking advantage of the constant component:
    simulate <- function (n, logpars, type, ub)
    {
        sigma <- exp(logpars[[1L]])
        d <- exp(logpars[[2L]])
        ## Sampling via polar coordinates and inversion method

        ## random angle
        theta <- runif(n, 0, 2*pi)

        ## sampling radius r
        ## trivial case u < sigma: p(r) \propto r on [0;u]
        if (ub < sigma) {
            r <- ub * sqrt(runif(n)) # inversion sampling
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
        unir <- runif(n) <= mass1 / (mass1 + mass2)

        ## samples from the uniform short-range component:
        n1 <- sum(unir)
        r1 <- sigma * sqrt(runif(n1)) # similar to the case u < sigma

        ## samples from power-law component: p2(r) \propto r^(-d+1) on [sigma;u]
        ## For d>2 only, we could use VGAM::rpareto(n,sigma,d-2), d=1 is trivial
        n2 <- n - n1
        r2 <- if (d==1) runif(n2, sigma, ub) else { # inversion sampling
            P2inv <- if (d == 2) { function (z) ub^z * sigma^(1-z) } else {
                function (z) (z*ub^(2-d) + (1-z)*sigma^(2-d))^(1/(2-d))
            }
            P2inv(runif(n2))
        }

        ## put samples from both components together
        r <- c(r1, r2)

        ## now rotate each point by a random angle to cover all directions
        r * cbind(cos(theta), sin(theta))
    }
    environment(simulate) <- getNamespace("stats")

    ## return the kernel specification
    list(f=f, F=F, Fcircle=Fcircle, deriv=deriv, Deriv=Deriv,
         simulate=simulate, npars=2L, validpars=validpars)
}


## integrate x*f(x) from 0 to R (vectorized)
intrfr.powerlawL <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    pl <- which(R > sigma)
    upper <- R
    upper[pl] <- sigma
    res <- upper^2 / 2                  # integral over x*constant part
    xplint <- if (d == 2) log(R[pl]/sigma) else (R[pl]^(2-d)-sigma^(2-d))/(2-d)
    res[pl] <- res[pl] + sigma^d * xplint
    res
}

## integrate x * (df(x)/dlogsigma) from 0 to R (vectorized)
intrfr.powerlawL.dlogsigma <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    pl <- which(R > sigma)
    res <- numeric(length(R))
    xplint <- if (d == 2) log(R[pl]/sigma) else (R[pl]^(2-d)-sigma^(2-d))/(2-d)
    res[pl] <- d * sigma^d * xplint
    res
}
## local({ # validation via numerical integration -> tests/testthat/test-siafs.R
##     p <- function (r, sigma, d)
##         r * siaf.powerlawL()$deriv(cbind(r,0), log(c(sigma,d)))[,1L]
##     Pnum <- function (r, sigma, d) sapply(r, function (.r) {
##         integrate(p, 0, .r, sigma=sigma, d=d, rel.tol=1e-8)$value
##     })
##     r <- c(1,2,5,10,20,50,100)
##     dev.null <- sapply(c(1,2,1.6), function(d) stopifnot(isTRUE(
##         all.equal(intrfr.powerlawL.dlogsigma(r, log(c(3, d))), Pnum(r, 3, d)))))
## })

## integrate x * (df(x)/dlogd) from 0 to R (vectorized)
intrfr.powerlawL.dlogd <- function (R, logpars, types = NULL)
{
    sigma <- exp(logpars[[1L]])
    d <- exp(logpars[[2L]])
    pl <- which(R > sigma)
    res <- numeric(length(R))
    res[pl] <- if (d == 2) -(sigma*log(R[pl]/sigma))^2 else
        (sigma^d * R[pl]^(2-d) * (d-2)*d*log(R[pl]/sigma) -
         d*(sigma^2 - R[pl]^(2-d)*sigma^d)) / (d-2)^2
    res
}
