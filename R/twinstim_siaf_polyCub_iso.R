################################################################################
### C-Level Cubature of "siaf" over Polygonal Domains using 'polyCub_iso'
###
### Copyright (C) 2017,2020 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################

### construct a call using either .polyCub.iso or its C-version
.call.polyCub.iso <- function (intrfr_name, engine = "C")
{
    if (engine == "C") {
        call("siaf_polyCub_iso", quote(polydomain$bdry), intrfr_name,
             quote(siafpars), quote(list(...)))
    } else {
        call(".polyCub.iso", quote(polydomain$bdry), as.name(intrfr_name),
             quote(siafpars), center = c(0,0), control = quote(list(...)))
    }
}

## construct siaf$F function
siaf_F_polyCub_iso <- function (intrfr_name, engine = "C")
{
    F <- function (polydomain, f, siafpars, type, ...) {}
    body(F) <- .call.polyCub.iso(intrfr_name, engine)
    environment(F) <- getNamespace("surveillance")
    return(F)
}

## construct siaf$Deriv function
siaf_Deriv_polyCub_iso <- function (intrfr_names, engine = "C")
{
    Deriv <- function (polydomain, deriv, siafpars, type, ...) {}
    res_names <- paste0("res", seq_along(intrfr_names))
    calls <- mapply(
        FUN = function (intrfr_name, res_name)
            call("<-", as.name(res_name), .call.polyCub.iso(intrfr_name, engine)),
        intrfr_name = intrfr_names, res_name = res_names,
        SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
    result <- as.call(c(as.name("c"), lapply(res_names, as.name)))
    body(Deriv) <- as.call(c(as.name("{"), calls, result))
    environment(Deriv) <- getNamespace("surveillance")
    return(Deriv)
}

## 'polys' is a list of polygons in the form of owin$bdry
## 'intrfr_name' identifies the function used in the integrand
## 'pars' is a vector of parameters for "intrfr"
siaf_polyCub_iso <- function (polys, intrfr_name, pars, control = list())
{
    ## default control arguments for polyCub_iso / Rdqags
    ## similar to args(stats::integrate)
    control <- modifyList(
        list(subdivisions = 100L, rel.tol = .Machine$double.eps^0.25,
             stop.on.error = TRUE),
        control)
    if (is.null(control[["abs.tol"]]))
        control$abs.tol <- control$rel.tol
    ## integrate over each polygon
    ints <- lapply(X = polys, FUN = siaf_polyCub1_iso,
                   intrfr_code = INTRFR_CODE[intrfr_name], pars = pars,
                   subdivisions = control$subdivisions,
                   rel.tol = control$rel.tol,
                   abs.tol = control$abs.tol,
                   stop.on.error = control$stop.on.error)
    sum(unlist(ints, recursive = FALSE, use.names = FALSE))
}

## 'xypoly' is a list(x, y) of vertex coordinates (open)
siaf_polyCub1_iso <- function (xypoly, intrfr_code, pars,
                               subdivisions = 100L,
                               rel.tol = .Machine$double.eps^0.25, abs.tol = rel.tol,
                               stop.on.error = TRUE)
{
    if (length(xypoly[["y"]]) != (L <- length(xypoly[["x"]])))
        stop("xypoly$x and xypoly$y must have equal length")
    .C("C_siaf_polyCub1_iso",
       as.double(xypoly$x), as.double(xypoly$y), as.integer(L),
       as.integer(intrfr_code), as.double(pars),
       as.integer(subdivisions), as.double(abs.tol), as.double(rel.tol),
       as.integer(stop.on.error),
       value = double(1L), abserr = double(1L), neval = integer(1L),
       PACKAGE = "surveillance")$value
}

## integer codes are used to select the corresponding C-routine,
## see ../src/twinstim_siaf_polyCub_iso.c
INTRFR_CODE <- c(
    "intrfr.powerlaw" = 10L,
    "intrfr.powerlaw.dlogsigma" = 11L,
    "intrfr.powerlaw.dlogd" = 12L,
    "intrfr.student" = 20L,
    "intrfr.student.dlogsigma" = 21L,
    "intrfr.student.dlogd" = 22L,
    "intrfr.powerlawL" = 30L,
    "intrfr.powerlawL.dlogsigma" = 31L,
    "intrfr.powerlawL.dlogd" = 32L,
    "intrfr.gaussian" = 40L,
    "intrfr.gaussian.dlogsigma" = 41L,
    "intrfr.exponential" = 50L,
    "intrfr.exponential.dlogsigma" = 51L
)
