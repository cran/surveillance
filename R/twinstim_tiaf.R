################################################################################
### Temporal interaction functions for twinstim's epidemic component.
### Specific implementations are in separate files (e.g.: exponential, step).
###
### Copyright (C) 2009-2014 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################



#####################
### "Constructor" ###
#####################

tiaf <- function (g, G, deriv, Deriv, npars, validpars = NULL)
{
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
    validpars <- if (!haspars || is.null(validpars))
        NULL else match.fun(validpars)
    list(g = g, G = G,
         deriv = deriv, Deriv = Deriv,
         npars = npars, validpars = validpars)
}



#################################
### Constant temporal interaction
#################################

tiaf.constant <- function ()
{
    res <- list(
        g = as.function(alist(t=, pars=, types=, rep.int(1, length(t))), envir = .GlobalEnv),
        G = as.function(alist(t=, pars=, types=, t), envir = .GlobalEnv),
        npars = 0L
    )
    attr(res, "constant") <- TRUE
    res
}

