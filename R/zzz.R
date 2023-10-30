#######################################
### Hook functions for package start-up
#######################################

.onLoad <- function (libname, pkgname)
{
    ## initialize options
    reset.surveillance.options()
}

.onAttach <- function (libname, pkgname)
{
    VERSION <- packageVersion(pkgname, lib.loc=libname)
    packageStartupMessage("This is ", pkgname, " ", VERSION, "; ",
                          "see ", sQuote(paste0("package?", pkgname)), " or\n",
                          "https://surveillance.R-Forge.R-project.org/",
                          " for an overview.")

    if (!interactive()) { # particularly for R CMD check
        ## skip long examples, unless:
        allExamples <- nzchar(Sys.getenv("_R_SURVEILLANCE_ALL_EXAMPLES_"))
        .Options$allExamples$value <- allExamples
    }
}



###########################
### Little helper functions
###########################


### determines multiplicities in a matrix (or data frame)
### and returns unique rows with appended column of counts
### using spatstat's multiplicity methods

countunique <- function (x) unique(cbind(x, COUNT = multiplicity(x)))
