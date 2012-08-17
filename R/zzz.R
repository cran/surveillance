#######################################
### Hook functions for package start-up
#######################################

## .onLoad <- function (libname, pkgname)
## {
##   #Load the CIdata thing
##   #data("CIdata", package=pkgname)

##   #Read the table of the hypgeom_2F1 function for parameters c(1/3,2/3) and
##   #5/3 -- atm this is computed for the values seq(0,10,by=0.01) and 11:100
##   #Load the pre-evaluated Hypergeometric function for computing Anscombe residuals
##   #data("hypGeomSmall",package=pkgname)

##   # <- those data sets should not be loaded into .GlobalEnv and be visible to the user!
##   # moved CIdata to (internal) sysdata
## }

.onAttach <- function (libname, pkgname)
{      
    ## Startup message
    vrs <- packageDescription(pkgname, lib.loc = libname, fields = "Version", drop = TRUE)
    packageStartupMessage("This is ", pkgname, " ", vrs, ". ",
                          "For overview type ",
                          sQuote(paste("?", pkgname, sep="")), ".")
    
    ## License limitation for package gpclib
    packageStartupMessage(paste(
        "\n\tNote: polygon geometry computations related to",
        "\n\tthe \"epidataCS\" class in ", pkgname, " depend on",
        "\n\tthe package gpclib, which has a restricted licence.\n",
        #"(Free for non-commercial use; commercial use prohibited.)",
        sep=""), appendLF = FALSE)
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


### Function which modifies a list _call_ according to another one similar to
### what the function utils::modifyList (by Deepayan Sarkar) does for list objects
## modifyListcall is used by update.twinstim

is.listcall <- function (x)
{
    is.call(x) &&
    as.character(x[[1]]) %in% c("list", "alist")
}

modifyListcall <- function (x, val)
{
    stopifnot(is.listcall(x), is.listcall(val))
    xnames <- names(x)[-1]
    for (v in names(val)[nzchar(names(val))]) {
        x[[v]] <-
            if (v %in% xnames && is.listcall(x[[v]]) && is.listcall(val[[v]]))
                modifyListcall(x[[v]], val[[v]]) else val[[v]]
    }
    x
}
