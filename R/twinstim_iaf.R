################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Spatial and temporal interaction functions for twinstim's epidemic component
###
### Copyright (C) 2009-2012 Sebastian Meyer
### $Revision: 459 $
### $Date: 2012-11-21 09:52:42 +0100 (Mi, 21. Nov 2012) $
################################################################################


### Function returning specification of constant spatial interaction/dispersal

## ## to avoid notes in R CMD check ("no visible binding for global variable ...")
## if (getRversion() >= "2.15.1") utils::globalVariables("r")

siaf.constant <- function () {
    r <- s <- n <- ub <- "just cheating on codetools::checkUsage"
    res <- list(
           f = as.function(alist(s=, pars=NULL, types=NULL, rep.int(1,nrow(s))),
                           envir = .GlobalEnv),
           Fcircle = as.function(alist(r=, pars=NULL, type=NULL, pi*r^2),
                                 envir = .GlobalEnv),
           ## simulation will be handled specially in simEpidataCS, this is only
           ## included here for completeness
           simulate = as.function(alist(n=, pars=NULL, type=NULL, ub=,
                                        surveillance:::runifdisc(n, ub)),
                                  envir = .GlobalEnv),
           npars = 0L
    )
    attr(res, "constant") <- TRUE
    res
}


### Generates a list specifying a Gaussian spatial interaction function
## nTypes: determines the number of parameters=(log-)standard deviations of the
##   Gaussian kernel. In a multitype epidemic, the different types may share the
##   same spatial interaction function (type-invariant), in which case nTypes=1.
##   Otherwise nTypes should equal the number of event types of the epidemic, in
##   which case every type has its own variance parameter.
## logsd: logical indicating if the gaussian kernel should be reparametrized
##   such that the log-standard deviation is the parameter in question. This
##   avoids constrained optimisation (L-BFGS-B) or the use of 'validpars'.
## density: logical. If TRUE, the isotropic Gaussian density (on R^2) will not
##   be scaled to have maximum value of 1 at the mean c(0,0).
## effRangeMult: determines the effective range for numerical integration in
##   terms of multiples of the parameter, i.e. with effRangeMult=6 numerical
##   integration only considers the 6-sigma area around the event instead of the
##   whole observation region W.
## validpars: If logsd = FALSE, you should either use
##   constrained optimisation (L-BFGS-B) or set 'validpars' to function (pars)
##   pars > 0. 

siaf.gaussian <- function (nTypes = 1, logsd = TRUE, density = FALSE,
                           F.adaptive = TRUE, effRangeMult = 6, validpars = NULL)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)

    f <- function (s, pars, types) {}       # coordinate matrix s, length(types) = 1 or nrow(s)
    F <- if (F.adaptive) {
        function(polydomain, f, pars, type, adapt=0.1) {}
    } else siaf.fallback.F
    Fcircle <- function (r, pars, type) {}  # single radius and type
    effRange <- function (pars) {}
    deriv <- function (s, pars, types) {}   # coordinate matrix s, length(types) = 1 or nrow(s)
    Deriv <- function (polydomain, deriv, pars, type, nGQ = 20L) {} # single "owin" and type
    simulate <- function (n, pars, type, ub) {} # n=size of the sample,
                                                # type=single type,
                                                # ub=upperbound (unused here)

    ## if there is only one type, we set the default type(s) argument to 1
    ## (it is actually unused inside the functions)
    if (nTypes == 1L) {
        formals(f)$types <- formals(Fcircle)$type <- formals(deriv)$types <-
            formals(Deriv)$type <- formals(simulate)$type <- 1L
    }
    
    # helper expressions
    tmp1 <- if (logsd) expression(sds <- exp(pars)) else expression(sds <- pars)
    tmp1.1 <- if (nTypes==1L) expression(sd <- sds) else expression(sd <- sds[type])
    tmp2 <- c(
            expression(
                sLengthSquared <- rowSums(s^2),
                L <- length(sLengthSquared)
                ),
            if (nTypes == 1L) expression(sdss <- sds) else expression(
                types <- rep(types, length.out = L),
                sdss <- sds[types]
                )
            )

    # spatial interaction function
    body(f) <- as.call(c(as.name("{"),
        tmp1, tmp2,
        expression(fvals <- exp(-sLengthSquared/2/sdss^2)),
        if (density) expression(fvals / (2*pi*sdss^2)) else expression(fvals)
    ))

    # numerically integrate f over a polygonal domain
    if (F.adaptive) {
        body(F) <- as.call(c(as.name("{"),
            tmp1, tmp1.1,
            expression(
                eps <- adapt * sd,
                intf <- polyCub.midpoint(polydomain, f, pars, type, eps=eps),
                intf
                )
        ))
    }
    
    # calculate the integral of f over a circular domain around 0
    body(Fcircle) <- as.call(c(as.name("{"),
        tmp1, tmp1.1,
        expression(val <- pchisq((r/sd)^2, 2)), # cf. Abramowitz&Stegun formula 26.3.24
        if (!density) expression(val <- val * 2*pi*sd^2),
        expression(val)
    ))

    # effective integration range of f as a function of sd
    if (isScalar(effRangeMult)) {
        body(effRange) <- as.call(c(as.name("{"),
            tmp1,
            substitute(effRangeMult*sds)
        ))
    } else effRange <- NULL

    # derivative of f wrt pars
    derivexpr <- if (logsd) { # derive f wrt psi=log(sd) !!
        if (density) {
            quote(deriv[cbind(1:L,colidx)] <- exp(-frac) / pi/sdss^2 * (frac-1))
        } else {
            quote(deriv[cbind(1:L,colidx)] <- exp(-frac) * 2*frac)
        }
    } else { # derive f wrt sd !!
        if (density) {
            quote(deriv[cbind(1:L,colidx)] <- exp(-frac) / pi/sdss^3 * (frac-1))
        } else {
            quote(deriv[cbind(1:L,colidx)] <- exp(-frac) * 2*frac/sdss)
        }
    }
    derivexpr <- do.call("substitute", args=list(expr=derivexpr,
                         env=list(colidx=if (nTypes==1L) 1L else quote(types))))
    body(deriv) <- as.call(c(as.name("{"),
        tmp1, tmp2,
        expression(
            deriv <- matrix(0, L, length(pars)),
            frac <- sLengthSquared/2/sdss^2
            ),
        derivexpr,
        expression(deriv)
    ))

    # integrate 'deriv' over a polygonal domain
    body(Deriv) <- as.call(c(as.name("{"),
        ## Determine a = argmax(abs(deriv(c(x,0))))
        if (density) {
            expression(a <- 0)          # maximum absolute value is at 0
        } else {
            c(tmp1, tmp1.1,
              expression(
                  xrange <- polydomain$xrange,           # polydomain is a "owin"
                  a <- min(max(abs(xrange)), sqrt(2)*sd), # maximum absolute value
                  if (sum(xrange) < 0) a <- -a # is more of the domain left of 0?
                  )
              )
        },
        if (nTypes == 1L) {
            expression(deriv.type <- function (s) deriv(s, pars, 1L)[,1L,drop=TRUE])
        } else { # d f(s|type_i) / d sigma_{type_j} is 0 for i != j
            expression(deriv.type <- function (s) deriv(s, pars, type)[,type,drop=TRUE])
        },
        expression(int <- polyCub.SV(polydomain$bdry, deriv.type, nGQ=nGQ, alpha=a)),
        if (nTypes == 1L) expression(int) else expression(
            res <- numeric(length(pars)), # zeros
            res[type] <- int,
            res
            )
    ))

    # sampler
    body(simulate) <- as.call(c(as.name("{"),
        tmp1, tmp1.1,
        expression(matrix(stats::rnorm(2*n, mean=0, sd=sd), nrow=n, ncol=2L))
    ))

    ## set function environments to the global environment
    environment(f) <- environment(F) <- environment(Fcircle) <-
        environment(deriv) <- environment(Deriv) <-
            environment(simulate) <- .GlobalEnv
    if (is.function(effRange)) environment(effRange) <- .GlobalEnv

    ## return the kernel specification
    list(f=f, F=F, Fcircle=Fcircle, effRange=effRange, deriv=deriv, Deriv=Deriv,
         simulate=simulate, npars=nTypes, validpars=validpars)
}



### Implementation of an isotropic power law kernel, specifically the
### Lomax distribution, i.e. a shifted Pareto distribution with domain [0;Inf)

## quantile function of the Lomax distribution
## we could also use VGAM::qlomax (but this would be slightly slower)
qlomax <- function (p, scale, shape) {
    scale * ((1-p)^(-1/shape) - 1)
}

## density=FALSE returns standardized Lomax kernel, i.e. f(x) = f_Lomax(x)/f(0),
## such that the kernel function starts at 1. f_Lomax(0) = alpha / sigma
siaf.lomax <- function (nTypes = 1, logpars = TRUE, density = FALSE,
                        effRangeProb = 0.999, validpars = NULL)
{
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
                polyCub.SV(polydomain$bdry, deriv1, paridx=paridx,
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
        exp <- sigma / (alpha-1)        # only valid if alpha > 1
        ##r <- VGAM::rlomax(n, scale=sigma, shape3.q=alpha)
        ## NO, r must be sampled from a density proportional to r*dlomax(r)!!!
        ##curve(x * alpha/sigma * (1+x/sigma)^-(alpha+1), 0, 200)

        ## We need the primitive of the density, i.e. int_0^r x*dlomax(x) dx
        cumd <- if (alpha == 1) {
            function (r) sigma * (log(1+r/sigma) - r/(r+sigma))
        } else {
            function (r)
                ifelse(is.infinite(r) & alpha > 1, exp,
                       exp * (1 - (1+r/sigma)^-alpha * (r*alpha/sigma+1)))
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
            function (q) cumd(q) / exp
        }
        
        ## However, there is no closed form expression for the quantile function
        ## of that distribution (inverse CDF), so we have to use uniroot
        ## However, uniroot needs a finite upper bound
        stopifnot(is.finite(ub))
        QF <- function(p) uniroot(function(q) CDF(q)-p, lower=0, upper=ub)$root

        ## now sample r as QF(U), where U ~ U(0,1)
        r <- sapply(runif(n), QF)
        
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

## Fcircle_Lomax <- function(r, logpars) {
##     logsigma <- logpars[[1]]  # used "[[" to drop names
##     logalpha <- logpars[[2]]
##     sigma <- exp(logsigma)
##     alpha <- exp(logalpha)
##     logfofr <- logalpha + alpha*logsigma - (alpha+1) * log(r+sigma)
##     fofr <- exp(logfofr)
##     ## method 1:
##     Intfinvsq <- function (z, sigma, alpha) {
##         ## antiderivative of ((f^-1)(z))^2
##         if (alpha == 1) {
##             sigma * (log(z) - 4*sqrt(sigma*z) + sigma*z)
##         } else {
##             (alpha*sigma^alpha)^(2/(alpha+1)) *
##                 z^((alpha-1)/(alpha+1))*(alpha+1)/(alpha-1) -
##                     2*alpha^(1/(alpha+1))*sigma^(1+alpha/(alpha+1)) *
##                         z^(alpha/(alpha+1))*(alpha+1)/alpha + sigma^2*z
##         }}
##     intfinvsq2 <- Intfinvsq(alpha/sigma, sigma, alpha) -
##         Intfinvsq(fofr, sigma, alpha)
##     ## method 2 (faster):
##     intfinvsq <- if (alpha == 1) {
##         sigma.rsigma <- sigma / (r+sigma)
##         sigma*(-2*log(sigma.rsigma) -3 +sigma.rsigma*(4-sigma.rsigma))
##     } else {
##         intfinvsq.fof0 <- 2*sigma / (alpha-1)
##         intfinvsq.fofr <- (sigma/(sigma+r))^alpha *
##             (r*(alpha+1)*(alpha*r+2*sigma)+2*sigma^2) /
##                 ((alpha-1)*(r+sigma))
##         intfinvsq.fof0 - intfinvsq.fofr
##     }
##     ## check if identical
##     if (!isTRUE(all.equal(intfinvsq2, intfinvsq))) {
##         browser("the two integrations of lomax^-1 differ")
##     }
##     base <- pi * r^2 * fofr   # volume of cylinder up to height f(r)
##     base + pi*intfinvsq
## }


### naive defaults for the siaf specification

## numerical integration of f over a polygonal domain (single "owin" and type)
siaf.fallback.F <- function(polydomain, f, pars, type, method = "SV", ...)
{
    if (identical(method,"SV"))
        polyCub.SV(polydomain, f, pars, type, alpha=0, ...) # since max at origin
    else 
        polyCub(polydomain, f, method, pars, type, ...)
}

## numerical integration of deriv over a polygonal domain
siaf.fallback.Deriv <- function (polydomain, deriv, pars, type, method = "SV", ...)
{
    deriv1 <- function (s, paridx)
        deriv(s, pars, type)[,paridx,drop=TRUE]
    intderiv1 <- function (paridx)
        polyCub(polydomain, deriv1, method, paridx=paridx, ...)
    derivInt <- sapply(seq_along(pars), intderiv1)
    derivInt
}



### Returns a list specifying an exponential temporal interaction function
# nTypes: cf. parameter description for siaf.gaussian

tiaf.exponential <- function (nTypes = 1)
{
    # time points vector t, length(types) = 1 or length(t)
    g <- function (t, alpha, types) {
        types <- rep(types, length.out = length(t))
        alphat <- alpha[types]
        exp(-alphat*t)
    }
    G <- function (t, alpha, types) {
        types <- rep(types, length.out = length(t))
        alphat <- alpha[types]
        ifelse(alphat==0, t, -exp(-alphat*t)/alphat)
    }
    deriv <- function (t, alpha, types) {
        L <- length(t)
        types <- rep(types, length.out = L)
        alphat <- alpha[types]
        deriv <- matrix(0, L, length(alpha))
        deriv[cbind(1:L,types)] <- -t*exp(-alphat*t)
        deriv
    }
    Deriv <- function (t, alpha, types) {
        L <- length(t)
        types <- rep(types, length.out = L)
        alphat <- alpha[types]
        Deriv <- matrix(0, L, length(alpha))
        Deriv[cbind(1:L,types)] <- ifelse(alphat==0, -t^2/2, (t+1/alphat)*exp(-alphat*t)/alphat)
        Deriv
    }

    list(g=g, G=G, deriv=deriv, Deriv=Deriv, npars=nTypes, validpars=NULL)
}


### Returns specification of constant temporal interaction

tiaf.constant <- function () {
    res <- list(
        g = as.function(alist(t=, pars=, types=, rep.int(1, length(t))), envir = .GlobalEnv),
        G = as.function(alist(t=, pars=, types=, t), envir = .GlobalEnv),
        deriv = NULL, Deriv = NULL, npars = 0L, validpars = NULL
    )
    attr(res, "constant") <- TRUE
    res
}



### Helper function which checks if the (spatial or temporal)
### interaction function 'iaf' returns 1 at the origin 0 or (0,0)

.checkiaf0 <- function(iaf, npars, validpars, type = c("siaf", "tiaf"))
{
    type <- match.arg(type)
    testpars <- rep.int(0.5, npars)
    if (!is.null(validpars)) {
        valid <- validpars(testpars)
        if (length(valid) != 1 || !is.logical(valid)) {
            stop("the '", type, "'/'validpars' function must return a single ",
                 "logical value")
        }
        if (!valid) {
            testpars <- rep.int(pi, npars)
            valid <- validpars(testpars)
        }
    } else valid <- TRUE
    if (valid) {
        ctext <- switch(type, siaf = "coordinate (0,0)", tiaf = "time point 0")
        iaf0 <- switch(type,
            siaf = iaf(t(c(0,0)), testpars, 1),
            tiaf = iaf(0, testpars, 1))
        if (!isTRUE(all.equal(iaf0, 1))) {
            message("CAVE: the '", type, "' function does not evaluate to 1 ",
                    "at ", ctext)
        }
    } else {
        message("the maximum value returned by the '", type,
                "' function should be 1")
    }
    NULL
}



### Checks if FUN has three arguments (s/t, pars, type) and
### eventually adds the last two

.checknargs3 <- function (FUN, name)
{
    FUN <- match.fun(FUN)
    NARGS <- length(formals(FUN))
    if (NARGS == 0L) {
        stop("the function '", name, "' must accept at least one argument")
    } else if (NARGS == 1L) {
        formals(FUN) <- c(formals(FUN), alist(pars=, types=))
    } else if (NARGS == 2L) {
        formals(FUN) <- c(formals(FUN), alist(types=))
    }
    FUN
}



### Checks the siaf specification
# f: spatial interaction function (siaf). must accept two arguments, a coordinate matrix and a parameter vector. for marked twinstim, it must accept a third argument which is the type of the event (either a single type for all locations or separate types for each location)
# F: function that integrates 'f' (2nd argument) over a polygonal domain (of class "owin", 1st argument). The third and fourth arguments are the parameters and the type, respectively. There may be additional arguments which are passed by the control.siaf list in twinstim().
# Fcircle: optional function for fast calculation of the integral of f over a circle with radius r (first argument). Further arguments like f. It must not be vectorized for model fitting (will be passed single radius and single type).
# effRange: optional function returning the effective range (domain) of f for the given set of parameters such that the circle with radius effRange contains the numerical essential proportion the integral mass, e.g. function (sigma) 6*sigma. The return value must be a vector of length nTypes (effRange for each type). Must be supplied together with Fcircle.
# deriv: optional derivative of f with respect to the parameters. Takes the same arguments as f but returns a matrix with as many rows as there were coordinates in the input and npars columns. The derivative is necessary for the score function.
# Deriv: function which integrates 'deriv' (2nd argument) over a polygonal domain (of class "owin", 1st argument). The third and fourth argument are the parameters and the type, respectively. There may be additional arguments which are passed by the control.siaf list in twinstim().
# simulate: optional function returning a sample drawn from the spatial kernel (i.e. a two-column matrix of 'n' points). The arguments are 'n' (size of the sample), 'pars' (parameter vector of the spatial kernel), for marked twinstim also 'type' (a single type of event being generated), and optionally 'ub' (upper bound, truncation of the kernel)
# npars: number of parameters
# validpars: a function indicating if a specific parameter vector is valid. Not necessary if npars == 0. If missing or NULL, it will be set to function (pars) TRUE. This function is rarely needed in practice, because usual box constrained parameters can be taken into account by using L-BFGS-B as the optimization method (with arguments 'lower' and 'upper').
# knots: not implemented. Knots (> 0) of a spatial interaction STEP function of the distance

checksiaf <- function (f, F, Fcircle, effRange, deriv, Deriv, simulate, npars, validpars, knots)
{
    # if siaf is a step function specified by knots
    if (!missing(knots)) {
        return(sort(unique(as.vector(knots,mode="numeric")), decreasing=FALSE))
    }

    # if siaf is a continuous function
    npars <- as.integer(npars)
    if (length(npars) != 1 || npars < 0L) {
        stop("'siaf$npars' must be a single nonnegative number")
    }
    f <- .checknargs3(f, "siaf$f")
    F <- if (missing(F) || is.null(F)) siaf.fallback.F else {
        F <- match.fun(F)
        if (length(formals(F)) < 4L)
            stop("siaf$F() must accept >=4 arguments ",
                 "(polydomain, f, pars, type)")
        F
    }
    haspars <- npars > 0L
    if (!haspars || missing(deriv)) deriv <- NULL
    if (!is.null(deriv)) deriv <- .checknargs3(deriv, "siaf$deriv")
    if (missing(effRange)) effRange <- NULL
    if (missing(Fcircle) || is.null(Fcircle)) {
        Fcircle <- NULL
        if (!is.null(effRange)) {
            message("'siaf$effRange' only works in conjunction with 'siaf$Fcircle'")
            effRange <- NULL
        }
    }
    if (!is.null(Fcircle)) Fcircle <- .checknargs3(Fcircle, "siaf$Fcircle")
    if (!is.null(effRange)) {
        effRange <- match.fun(effRange)
        if (length(formals(effRange)) < 1L) {
            stop("the 'siaf$effRange' function must accept a parameter vector")
        }
    }
    Deriv <- if (is.null(deriv)) NULL else if (missing(Deriv) || is.null(Deriv))
        siaf.fallback.Deriv else {
            Deriv <- match.fun(Deriv)
            if (length(formals(Deriv)) < 4L)
                stop("siaf$Deriv() must accept >=4 arguments ",
                     "(polydomain, deriv, pars, type)")
            Deriv
        }
    ## Check if simulation function has proper format
    if (missing(simulate)) simulate <- NULL
    if (!is.null(simulate)) {
        simulate <- .checknargs3(simulate, "siaf$simulate")
        if (length(formals(simulate)) == 3L)
            formals(simulate) <- c(formals(simulate), alist(ub=))
    }
    ## Check if the validpars are of correct form
    validpars <- if (!haspars || missing(validpars) || is.null(validpars))
        NULL else match.fun(validpars)
    ## Check if the siaf has value 1 at (0,0)
    .checkiaf0(f, npars, validpars, "siaf")
    ## Done, return result.
    list(f = f, F = F, Fcircle = Fcircle, effRange = effRange,
         deriv = deriv, Deriv = Deriv,
         simulate = simulate,
         npars = npars, validpars = validpars)
}



### Checks the tiaf specification
# g: temporal interaction function (tiaf). must accept two arguments, a vector of time points and a parameter vector. for marked twinstim, it must accept a third argument which is the type of the event (either a single type for all locations or separate types for each location)
# G: a primitive of g. Must accept same arguments as g. (_vector_ of time points, not just a single one!)
# deriv: optional derivative of g with respect to the parameters. Takes the same arguments as g but returns a matrix with as many rows as there were time points in the input and npars columns. The derivative is necessary for the score function.
# Deriv: optional primitive of deriv (with respect to time). Must accept same arguments as deriv and g. Returns a matrix with as many rows as there were time points in the input and npars columns. The integrated derivative is necessary for the score function.
# npars: number of parameters, which are passed as second argument to g and G
# validpars: a function indicating if a specific parameter vector is valid. Not necessary if npars == 0. If missing or NULL, will be set to function (pars) TRUE. This function is rarely needed in practice, because usual box constrained parameters can be taken into account by using L-BFGS-B as the optimization method (with arguments 'lower' and 'upper').
# knots: not implemented. Knots (> 0) of a temporal interaction STEP function

checktiaf <- function (g, G, deriv, Deriv, npars, validpars, knots)
{
    # if tiaf is a step function specified by knots
    if (!missing(knots)) {
        stop("'knots' are not implemented for 'tiaf'")
        return(sort(unique(as.vector(knots,mode="numeric")), decreasing=FALSE))
    }

    # if tiaf is a continous function
    npars <- as.integer(npars)
    if (length(npars) != 1 || npars < 0L) {
        stop("'tiaf'/'npars' must be a single nonnegative number")
    }
    haspars <- npars > 0L
    g <- .checknargs3(g, "tiaf$g")
    G <- .checknargs3(G, "tiaf$G")
    if (!haspars || missing(deriv)) deriv <- NULL
    if (!haspars || missing(Deriv)) Deriv <- NULL
    if (!is.null(deriv)) deriv <- .checknargs3(deriv, "tiaf$deriv")
    if (!is.null(Deriv)) Deriv <- .checknargs3(Deriv, "tiaf$Deriv")
    validpars <- if (!haspars || missing(validpars) || is.null(validpars))
        NULL else match.fun(validpars)
    ## Check if the tiaf has value 1 at the origin
    .checkiaf0(g, npars, validpars, "tiaf")
    list(g = g, G = G, deriv = deriv, Deriv = Deriv, npars = npars, validpars = validpars)
}



### Wrapper used in main function to evaluate the siaf and tiaf arguments.
### If succesful, returns checked interaction function

.parseiaf <- function (iaf)
{
    name <- as.character(substitute(iaf))
    if (missing(iaf) || is.null(iaf)) {
        message("assuming constant ",
                switch(name, siaf="spatial", tiaf="temporal"),
                " interaction '", name, ".constant()'")
        do.call(paste(name, "constant", sep="."), args=alist())
    } else if (is.list(iaf)) {
        ret <- do.call(paste0("check",name), args = iaf)
        attr(ret, "constant") <- isTRUE(attr(iaf, "constant"))
        ret
    } else if (is.vector(iaf, mode = "numeric")) {
        stop("'knots' are not implemented for '",name,"'")
        do.call(paste0("check",name), args = list(knots = iaf))
    } else {
        stop("'",name,"' must be NULL (or missing), a list (-> continuous ",
            "function), or numeric (-> knots of step function)")
    }
}
