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
    allExamples <- if (interactive()) {
        TRUE
    } else { # R CMD check
        ## only do all examples if a specific environment variable is set
        ## (to any value different from "")
        nzchar(Sys.getenv("_R_SURVEILLANCE_ALL_EXAMPLES_"))
        ## CAVE: testing for _R_CHECK_TIMINGS_ as in surveillance < 1.9-1
        ## won't necessarily skip long examples for daily checks on CRAN (see
        ## https://stat.ethz.ch/pipermail/r-devel/2012-September/064812.html
        ## ). For instance, the daily Windows checks run without timings.
    }
    surveillance.options(allExamples = allExamples)
}



###########################
### Little helper functions
###########################


### determines multiplicities in a matrix (or data frame)
### and returns unique rows with appended column of counts
### using spatstat's multiplicity methods

countunique <- function (x) unique(cbind(x, COUNT = multiplicity(x)))
