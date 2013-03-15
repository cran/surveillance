#######################################
### Hook functions for package start-up
#######################################

### not sure if we should query gpclib permission via spatstat or maptools:
### within a single R session, spatstat might be used for commercial purposes
### _without_ gpclib, whereas some functionality of the surveillance package might
### be used non-commercially _with_ gpclib
## gpclibPermitStatus <- function ()
## {
##     ## check global gpclib permission status
##     ##globally <- isTRUE(getOption("gpclib"))
##     ## -> CRAN team does not like packages which use own global options
##
##     ## check for permission via surveillance.options
##     via.surveillance <- surveillance.options("gpclib")
##
##     ## check for permission via spatstat package
##     via.spatstat <- spatstat.options("gpclib")
##
##     ## check for permission via maptools
##     via.maptools <- if ("maptools" %in% loadedNamespaces())
##         maptools::gpclibPermitStatus() else FALSE
##
##     ## return gpclib permission status
##     via.surveillance | via.spatstat | via.maptools
## }

gpclibCheck <- function (fatal = TRUE)
{
    gpclibOK <- surveillance.options("gpclib")
    if (!gpclibOK && fatal) {
        message("Note: The gpclib license is accepted by ",
                sQuote("surveillance.options(gpclib=TRUE)"), ".")
        stop("acceptance of the gpclib license is required")
    }
    gpclibOK
}

.onLoad <- function (libname, pkgname)
{
    ## Determine package version and store it as surveillance:::VERSION
    v <- packageDescription(pkgname, lib.loc=libname, fields="Version", drop=TRUE)
    ##<- a convenience function packageVersion() was only introduced in R 2.12.0
    assign("VERSION", package_version(v), getNamespace(pkgname))

    ## initialize options
    reset.surveillance.options()
}

.onAttach <- function (libname, pkgname)
{
    ## Startup message
    packageStartupMessage("This is ", pkgname, " ", VERSION, ". ",
                          "For overview type ",
                          sQuote(paste0("help(", pkgname, ")")), ".")

    ## License limitation for package gpclib
    packageStartupMessage(
        "Note: Polygon geometry computations required for the generation of",
        "\n      \"epidataCS\" objects currently depend on the gpclib package,",
        "\n      which has a restricted license.  This functionality is disabled",
        "\n      by default, but may be enabled by setting",
        "\n      ", sQuote("surveillance.options(gpclib=TRUE)"),
                 ", if this license is applicable."
    )

    ## decide if we should run all examples (some take a few seconds)
    allExamples <- if (interactive()) TRUE else {
        .withTimings <- Sys.getenv("_R_CHECK_TIMINGS_")
        withTimings <- !is.na(.withTimings) && nzchar(.withTimings)
        !withTimings
    }
    surveillance.options(allExamples = allExamples)
}


### Function 'base::paste0()' only exists as of R version 2.15.0
### Define it as a wrapper for base::paste() for older versions

if (getRversion() < "2.15.0" || R.version$"svn rev" < 57795 ||
    !exists("paste0", baseenv())) {
    paste0 <- function (..., collapse = NULL) {
        ## the naive way: paste(..., sep = "", collapse = collapse)
        ## probably better: establish appropriate paste() call:
        cl <- match.call()
        names(cl) <- sub("sep", "", names(cl)) # no sep argument
        cl$sep <- ""
        cl[[1]] <- as.name("paste")
        eval(cl, envir = parent.frame())
    }
}



###########################
### Little helper functions
###########################


### checking if x is scalar, i.e. a numeric vector of length 1.

isScalar <- function (x) {
    length(x) == 1L && is.vector(x, mode = "numeric")
}


### returns the dot/scalar product of two vectors

dotprod <- function (x,y)
{
    sum(x*y)
}


### _c_onditional lapply, which only uses lapply() if X really is a list object
### and otherwise applies FUN to X. The result is always a list (of length 1 in
### the latter case). Used for neOffset in hhh4 models.

clapply <- function (X, FUN, ...)
{
    if (is.list(X)) lapply(X, FUN, ...) else list(FUN(X, ...))
}


### pretty p-value formatting

formatPval <- function (pv, eps = 1e-4)
{
    format1 <- function (p)
        format.pval(p, digits = if (p<10*eps) 1 else 2, eps = eps)
    sapply(pv, format1)
}


### quantile function of the Lomax distribution
### we could also use VGAM::qlomax (but this would be slightly slower)

qlomax <- function (p, scale, shape) {
    scale * ((1-p)^(-1/shape) - 1)
}
