################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Temporal interaction functions for twinstim's epidemic component.
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 527 $
### $Date: 2013-03-08 13:25:48 +0100 (Fr, 08. Mrz 2013) $
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
