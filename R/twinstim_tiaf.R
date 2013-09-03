################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Temporal interaction functions for twinstim's epidemic component.
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 542 $
### $Date: 2013-04-28 17:03:50 +0200 (Son, 28 Apr 2013) $
################################################################################



#####################
### "Constructor" ###
#####################

## g: temporal interaction function (tiaf). must accept two arguments, a vector
##    of time points and a parameter vector. for marked twinstim, it must accept
##    a third argument which is the type of the event (either a single type for
##    all locations or separate types for each location).
## G: a primitive of g. Must accept same arguments as g. (_vector_ of time
##    points, not just a single one!) 
## deriv: optional derivative of g with respect to the parameters. Takes the
##    same arguments as g but returns a matrix with as many rows as there were
##    time points in the input and npars columns. The derivative is necessary
##    for the score function. 
## Deriv: optional primitive of deriv (with respect to time). Must accept same
##    arguments as deriv and g. Returns a matrix with as many rows as there were
##    time points in the input and npars columns. The integrated derivative is
##    necessary for the score function. 
## npars: number of parameters, which are passed as second argument to g and G
## validpars: optional function indicating if a specific parameter vector is
##    valid. If missing or NULL, will be set to function (pars) TRUE.
##    This function is rarely needed in practice, because usual box
##    constrained parameters can be taken into account by using L-BFGS-B as the
##    optimization method (with arguments 'lower' and 'upper'). 
## knots: not implemented. Knots (> 0) of a temporal interaction STEP function

tiaf <- function (g, G, deriv, Deriv, npars, validpars, knots)
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
    list(g = g, G = G, deriv = deriv, Deriv = Deriv, npars = npars, validpars = validpars)
}



#################################
### Constant temporal interaction
#################################

tiaf.constant <- function ()
{
    res <- list(
        g = as.function(alist(t=, pars=, types=, rep.int(1, length(t))), envir = .GlobalEnv),
        G = as.function(alist(t=, pars=, types=, t), envir = .GlobalEnv),
        deriv = NULL, Deriv = NULL, npars = 0L, validpars = NULL
    )
    attr(res, "constant") <- TRUE
    res
}



#############################################
### Exponential temporal interaction function
#############################################

## nTypes: determines the number of parameters=(log-)alphas of the
##    Exponential kernel. In a multitype epidemic, the different types may share
##    the same temporal interaction function (type-invariant), in which case
##    nTypes=1. Otherwise nTypes should equal the number of event types of the
##    epidemic, in which case every type has its own alpha.

tiaf.exponential <- function (nTypes = 1)
{
    nTypes <- as.integer(nTypes)
    stopifnot(length(nTypes) == 1L, nTypes > 0L)
    
    ## function definitions for nTypes = 1 (length(alpha) == 1)
    g <- function (t, alpha, types) {
        exp(-alpha*t)
    }
    G <- function (t, alpha, types) {
        if (alpha==0) t else -exp(-alpha*t)/alpha
    }
    deriv <- function (t, alpha, types) {
        as.matrix( -t*exp(-alpha*t) )
    }
    Deriv <- function (t, alpha, types) {
        as.matrix( if (alpha==0) -t^2/2 else (t+1/alpha)*exp(-alpha*t)/alpha )
    }

    ## adaptions for nTypes > 1
    if (nTypes > 1) {
        ## time points vector t, length(types) = length(t)
        body(g) <- as.call(append(as.list(body(g)),
                                  quote(alpha <- alpha[types]), after=1))
        body(G) <- quote({
            alpha <- alpha[types]
            ifelse (alpha==0, t, -exp(-alpha*t)/alpha)
        })
        body(deriv) <- quote({
            L <- length(t)
            deriv <- matrix(0, L, length(alpha))
            alpha <- alpha[types]
            deriv[cbind(1:L,types)] <- -t*exp(-alpha*t)
            deriv
        })
        body(Deriv) <- quote({
            L <- length(t)
            Deriv <- matrix(0, L, length(alpha))
            alpha <- alpha[types]
            Deriv[cbind(1:L,types)] <-
                ifelse(alpha==0, -t^2/2,
                       (t+1/alpha)*exp(-alpha*t)/alpha)
            Deriv
        })
    }

    ## set function environments to the global environment
    environment(g) <- environment(G) <-
        environment(deriv) <- environment(Deriv) <- .GlobalEnv

    ## return the kernel specification
    list(g=g, G=G, deriv=deriv, Deriv=Deriv, npars=nTypes, validpars=NULL)
}
