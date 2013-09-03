################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Spatial interaction functions for twinstim's epidemic component.
### Specific implementations are in seperate files (e.g.: Gaussian, power law).
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 627 $
### $Date: 2013-08-22 22:39:06 +0200 (Don, 22 Aug 2013) $
################################################################################



#####################
### "Constructor" ###
#####################

## f: spatial interaction function (siaf). must accept two arguments, a
##    coordinate matrix and a parameter vector. for marked twinstim, it must
##    accept a third argument which is the type of the event (either a single
##    type for all locations or separate types for each location).
## F: function that integrates 'f' (2nd argument) over a polygonal domain (of
##    class "owin", 1st argument). The third and fourth arguments are the
##    parameters and the type, respectively. There may be additional arguments
##    which are passed by the control.siaf list in twinstim(). 
## Fcircle: optional function for fast calculation of the integral of f over a
##    circle with radius r (first argument). Further arguments like f. It
##    must not be vectorized for model fitting (will be passed single
##    radius and single type). 
## effRange: optional function returning the effective range (domain) of f for
##    the given set of parameters such that the circle with radius effRange
##    contains the numerical essential proportion the integral mass, e.g.
##    function (sigma) 6*sigma. The return value must be a vector of length
##    nTypes (effRange for each type). Must be supplied together with Fcircle. 
## deriv: optional derivative of f with respect to the parameters. Takes the
##    same arguments as f but returns a matrix with as many rows as there were
##    coordinates in the input and npars columns. The derivative is necessary
##    for the score function. 
## Deriv: function which integrates 'deriv' (2nd argument) over a polygonal
##    domain (of class "owin", 1st argument). The third and fourth argument are
##    the parameters and the type, respectively. There may be additional
##    arguments which are passed by the control.siaf list in twinstim(). 
## simulate: optional function returning a sample drawn from the spatial kernel
##    (i.e. a two-column matrix of 'n' points). The arguments are 'n' (size of
##    the sample), 'pars' (parameter vector of the spatial kernel), for marked
##    twinstim also 'type' (a single type of event being generated), and
##    optionally 'ub' (upper bound, truncation of the kernel).
## npars: number of parameters
## validpars: optional function indicating if a specific parameter vector is
##    valid. If missing or NULL, it will be set to function
##    (pars) TRUE. This function is rarely needed in practice, because usual box
##    constrained parameters can be taken into account by using L-BFGS-B as the
##    optimization method (with arguments 'lower' and 'upper'). 
## knots: not implemented. Knots (> 0) of a spatial interaction STEP function of the distance

siaf <- function (f, F, Fcircle, effRange, deriv, Deriv, simulate, npars, validpars, knots)
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
    ## Done, return result.
    list(f = f, F = F, Fcircle = Fcircle, effRange = effRange,
         deriv = deriv, Deriv = Deriv,
         simulate = simulate,
         npars = npars, validpars = validpars)
}



##########################################
### Constant spatial interaction/dispersal
##########################################

siaf.constant <- function ()
{
    r <- s <- n <- ub <- "just cheating on codetools::checkUsage"
    ## to avoid notes in R CMD check ("no visible binding for global variable")
    ## one could also use utils::globalVariables() in R >= 2.15.1 as follows:
    ## if (getRversion() >= "2.15.1") utils::globalVariables("r")
    res <- list(
           f = as.function(alist(s=, pars=NULL, types=NULL, rep.int(1,nrow(s))),
                           envir = .GlobalEnv),
           Fcircle = as.function(alist(r=, pars=NULL, type=NULL, pi*r^2),
                                 envir = .GlobalEnv),
           ## simulation will be handled specially in simEpidataCS, this is only
           ## included here for completeness
           simulate = as.function(alist(n=, pars=NULL, type=NULL, ub=,
                                        runifdisc(n, ub)),
                                  envir = getNamespace("surveillance")),
           npars = 0L
    )
    attr(res, "constant") <- TRUE
    res
}



##########################################
### naive defaults for the siaf primitives
##########################################

## numerical integration of f over a polygonal domain (single "owin" and type)
siaf.fallback.F <- function(polydomain, f, pars, type, method = "SV", ...)
{
    if (identical(method,"SV"))
        polyCub::polyCub.SV(polydomain, f, pars, type, alpha=0, ...) # since max at origin
    else 
        polyCub::polyCub(polydomain, f, method, pars, type, ...)
}

## numerical integration of deriv over a polygonal domain
siaf.fallback.Deriv <- function (polydomain, deriv, pars, type, method = "SV", ...)
{
    deriv1 <- function (s, paridx)
        deriv(s, pars, type)[,paridx,drop=TRUE]
    intderiv1 <- function (paridx)
        polyCub::polyCub(polydomain, deriv1, method, paridx=paridx, ...)
    derivInt <- sapply(seq_along(pars), intderiv1)
    derivInt
}



######################################
### Check Fcircle, deriv, and simulate
######################################

checksiaf <- function (siaf, pargrid, type=1)
{
    stopifnot(is.list(siaf), is.numeric(pargrid), !is.na(pargrid),
              length(pargrid) > 0)
    pargrid <- as.matrix(pargrid)
    stopifnot(siaf$npars == ncol(pargrid))
    
    ## Check 'Fcircle'
    if (!is.null(siaf$Fcircle)) {
        cat("'Fcircle' vs. numerical cubature ... ")
        res <- checksiaf.Fcircle(siaf$Fcircle, siaf$f, pargrid, type=type)
        cat(all.equal(res[,1],res[,2]), "\n")
    }
    
    ## Check 'deriv'
    if (!is.null(siaf$deriv)) {
        cat("'deriv' vs. numerical derivative ... ")
        maxRelDiffs <- checksiaf.deriv(siaf$deriv, siaf$f, pargrid, type=type)
        cat("maxRelDiff =", max(maxRelDiffs), "\n")
    }

    ## Check 'simulate'
    if (!is.null(siaf$simulate)) {
        checksiaf.simulate(siaf$simulate, siaf$f, pargrid[1,], type=type)
    }
}

checksiaf.Fcircle <- function (Fcircle, f, pargrid, type=1,
                               rs=runif(20, 0, 100), nGQ=30)
{
    pargrid <- pargrid[rep(1:nrow(pargrid), each=length(rs)),,drop=FALSE]
    rpargrid <- cbind(rs, pargrid, deparse.level=0)
    res <- t(apply(rpargrid, 1, function (x) {
        c(ana = Fcircle(x[1], x[-1], type),
          num = polyCub.SV(discpoly(c(0,0), x[1], npoly=128, class="owin"),
                           function (s) f(s, x[-1], type),
                           alpha=0, nGQ=nGQ))
    }))
    res
}

checksiaf.deriv <- function (deriv, f, pargrid, type=1, rmax=100)
{
    rgrid <- seq(-rmax,rmax,len=21) / sqrt(2)
    rgrid <- rgrid[rgrid != 0] # some siafs are always 1 at (0,0) (deriv=0)
    sgrid <- cbind(rgrid, rgrid)
    maxreldiffs <- if (requireNamespace("maxLik")) {
        apply(pargrid, 1, function (pars) {
            maxLik::compareDerivatives(f, deriv, t0=pars, s=sgrid,
                                       print=FALSE)$maxRelDiffGrad
        })
    } else { # use stats::numericDeriv
        apply(pargrid, 1, function (pars) {
            ana <- deriv(sgrid, pars, types=type)
            logsigma <- pars[1L]; logd <- pars[2L]
            num <- attr(stats::numericDeriv(quote(f(sgrid, c(logsigma,logd),
                                                    types=type)),
                                            theta=c("logsigma", "logd")),
                        "gradient")
            max((ana-num)/(0.5*(abs(ana)+abs(num))))
        })
    }
    maxreldiffs
}

checksiaf.simulate <- function (simulate, f, pars, type=1, B=3000, ub=10)
{
    ## Simulate B points on the disc with radius 'ub'
    simpoints <- simulate(B, pars, type=type, ub=ub)

    ## Graphical check
    par(mar=c(1,2,2,1))
    plot(spatstat::as.im.function(function(x,y,...) f(cbind(x,y), pars, type),
                                  W=discpoly(c(0,0), ub, class="owin")),
         axes=TRUE, main="Simulation from the spatial kernel")
    points(simpoints, cex=0.2)
    kdens <- MASS::kde2d(simpoints[,1], simpoints[,2], n=100)
    contour(kdens, add=TRUE, col=2, lwd=2,
            labcex=1.5, vfont=c("sans serif", "bold"))
    ##x11(); image(kdens, add=TRUE)
}
