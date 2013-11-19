#######################################
### Hook functions for package start-up
#######################################

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
    ## initialize options
    reset.surveillance.options()
}

.onAttach <- function (libname, pkgname)
{
    ## Startup message
    VERSION <- packageVersion(pkgname, lib.loc=libname)
    packageStartupMessage("This is ", pkgname, " ", VERSION, ". ",
                          "For overview type ",
                          sQuote(paste0("help(", pkgname, ")")), ".")

    ## decide if we should run all examples (some take a few seconds)
    allExamples <- if (interactive()) TRUE else {
        .withTimings <- Sys.getenv("_R_CHECK_TIMINGS_")
        withTimings <- !is.na(.withTimings) && nzchar(.withTimings)
        !withTimings
    }
    surveillance.options(allExamples = allExamples)
}


### Function 'base::rep_len' only exists as of R version 3.0.0
### Define it as a wrapper for base::rep() for older versions

if (getRversion() < "3.0.0" || !exists("rep_len", baseenv())) {
    rep_len <- function (x, length.out) rep(x, length.out=length.out)
}



###########################
### Little helper functions
###########################


### checking if x is scalar, i.e. a numeric vector of length 1.

isScalar <- function (x) {
    length(x) == 1L && is.vector(x, mode = "numeric")
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
