################################################################################
### Set up .Options for "surveillance".
### Inspired by the options management in the R package "spatstat" (1.29-0)
### by Adrian Baddeley and Rolf Turner (2012).
###
### Copyright (C) 2012,2023 Sebastian Meyer, 2014 Michael Hoehle
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################

.Options <- new.env()

## the "gpclib" option is obsolete and no longer documented
.Options$gpclib <- list(
    default = FALSE, # no gpclib due to license restrictions
    check = function(x) {
        if (!is.logical(x) || length(x) != 1L) return(FALSE)
        if (x) {
            .Deprecated(msg = paste("Option", dQuote("gpclib"), "is obsolete."))
            if (!requireNamespace("gpclib")) {
                warning("cannot use ", sQuote("gpclib"), call. = FALSE)
                return(FALSE)
            }
        }
        TRUE
    },
    valid = "a single logical value"
    )

## mainly a maintainer-level switch driven by _R_SURVEILLANCE_ALL_EXAMPLES_
.Options$allExamples <- list(
    default = TRUE,  # maybe disabled by .onAttach()
    check = function(x) is.logical(x) && length(x) == 1L,
    valid = "a single logical value"
    )

#Tick sizes of sts xaxis relative to par()$tcl
.Options$stsTickFactors <- list(
    default = c("%d"=0.33,"%W"=0.33,"%V"=0.33,"%m"=1,"%Q"=1.25,"%Y"=1.5,"%G"=1.5),
    check = function(x) is.vector(x, mode="numeric") && !is.null(names(x)),
    valid = "a named vector of relative tick sizes"
)

#Colors for the prediction intervals in nowcast plots
.Options$colors <- list(
    default = c(nowSymbol="springgreen4",piBars="orange"),
    check = function(x) is.character(x),
    valid = "a vector of color names"
)


## Function to activate the defaults
reset.surveillance.options <- function ()
{
    opts <- sapply(ls(.Options, all.names=TRUE), function (option) {
        .Options[[option]]$value <- .Options[[option]]$default
    }, simplify=FALSE, USE.NAMES=TRUE)
    invisible(opts)
}

## Internal function to query options
get.surveillance.options <- function (x, drop = TRUE)
{
    opts <- lapply(.Options, "[[", "value")
    if (drop && !missing(x) && length(x) == 1L) opts[[x]] else opts[x]
}

## Exported function to modify and query options
surveillance.options <- function (...) 
{
    knownOptions <- ls(.Options, all.names=TRUE)
    
    called <- list(...)
    if (length(called) == 0) return(get.surveillance.options())
    if (is.null(names(called)) && length(called)==1) {
      x <- called[[1]]
      if (is.null(x)) return(get.surveillance.options())
      if (is.list(x)) called <- x
    }
    
    if (is.null(names(called))) # case: surveillance.options("par1","par2",...)
    {
	ischar <- unlist(lapply(called, is.character))
	if(all(ischar)) {
          choices <- unlist(called)
          ok <- choices %in% knownOptions
          if(!all(ok)) stop("unrecognised option(s): ", called[!ok])
          return(get.surveillance.options(choices))
	} else {
	   wrong <- called[!ischar]
	   offending <- unlist(lapply(wrong, deparse, nlines=1,
                                      control="delayPromises"))
	   offending <- paste(offending, collapse=",")
           stop("unrecognised mode of argument(s) [", offending, "]:",
                "\n  should be character string or name=value pair")
    	}
    } else { # case: surveillance.options(name=value, name2=value2, ...)
        assignto <- names(called)
        if (!all(nzchar(assignto)))
            stop("options must all be identified by name=value")
        recog <- assignto %in% knownOptions
        if(!all(recog))
            stop("unrecognised option(s): ", assignto[!recog])
        
        ## validate and assign new values
        oldopts <- get.surveillance.options(assignto, drop=FALSE)
        for(i in seq_along(assignto)) {
            nama <- assignto[i]
            valo <- called[[i]]
            entry <- .Options[[nama]]
            if (!entry$check(valo))
                stop("option ", dQuote(nama), " should be ", entry$valid)
            .Options[[nama]]$value <- valo
        }
        
        ## done
        invisible(oldopts)
    }
}
