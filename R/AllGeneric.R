
### Define some functions to be S3 generic

animate <- function (object, ...) UseMethod("animate")
R0 <- function (object, ...) UseMethod("R0")
as.epidata <- function (data, ...) UseMethod("as.epidata")
intensityplot <- function (x, ...) UseMethod("intensityplot")
untie <- function (x, amount, ...) UseMethod("untie")
intersectPolyCircle <- function (object, center, radius, ...)
    UseMethod("intersectPolyCircle")
calibrationTest <- function (x, ...) UseMethod("calibrationTest")
scores <- function (x, ...) {
    if (identical(class(x), "list")) {
        ## backward compatibility with surveillance < 1.10-0
        scores.oneStepAhead(x, ...)
    } else {
        UseMethod("scores")
    }
}
pit <- function (x, ...) UseMethod("pit")

## internal function with methods for "twinSIR" and "simEpidata"
getModel <- function (object, ...) UseMethod("getModel")

## list coefficients by component
coeflist <- function (x, ...) UseMethod("coeflist")
coeflist.default <- function (x, npars, ...)
{
    if (is.null(groupnames <- names(npars))) {
        stop("'npars' must be named")
    }
    f <- factor(rep.int(groupnames, npars), levels = groupnames)
    split.default(x = x, f = f, drop = FALSE)
}


### Declare some existing R functions (which we import) to be S4-generic.
### This is not strictly necessary, but considered better programming style, and
### it avoids messages noting the creation of the generics during package build
### and installation, see the section "Basic Use" in help("setGeneric").

setGeneric("plot")
setGeneric("aggregate")
setGeneric("toLatex")

## data frame-like methods defined in sts.R
setGeneric("dim")
setGeneric("dimnames")

## "ts"-like methods
setGeneric("frequency")
setGeneric("start")


######################################################################
# Conversion to and from sts objects
######################################################################

#setGeneric("as.sts")
setGeneric("as.data.frame")


######################################################################
# Accessing and replacing slots of the "sts" class
######################################################################
#epoch slot
setGeneric("epoch", function(x, as.Date=x@epochAsDate) standardGeneric("epoch"))
setMethod("epoch", "sts", function(x, as.Date=x@epochAsDate) {
  if (!as.Date) { # return numeric vector
    x@epoch
  } else { # convert to Date format
    if (x@epochAsDate) {
      as.Date(x@epoch, origin = "1970-01-01")
    } else if (x@freq == 12) { # use the first day of every month
      as.Date(strptime(paste(year(x), epochInYear(x), 1, sep = "-"),
                       format = "%Y-%m-%d"))
    } else if (x@freq %in% c(26, 52)) { # use Mondays
      ## be consistent with epochInYear(): 'start' means *ISO* year and week!
      ## Unfortunately, %G and %V are not supported for input via strptime():
      ## firstMonday <- strptime(x = paste0(x@start[1L], "-W", x@start[2L], "-1"),
      ##                         format = "%G-W%V-%u")  # WRONG, gives today
      ## so we run a naive search for the Monday of the 'start' week
      startweek <- x@start[2L]
      if (x@freq == 26) startweek <- 2*startweek - 1
      candidates <- seq.Date(as.Date(paste0(x@start[1L]-1L, "-12-29")),
                             as.Date(paste0(x@start[1L], "-12-28")), by = 1L)
      firstMonday <- candidates[match(sprintf("%i-W%02i", x@start[1L], startweek),
                                      strftime(candidates, "%G-W%V"))]
      seq.Date(from = firstMonday, by = (52/x@freq) * 7L, length.out = nrow(x))
    } else if (x@freq == 365) { # use day of the year (incorrect in leap years)
      as.Date(strptime(paste0(year(x), "-D", epochInYear(x)), format = "%Y-D%j"))
    } else {
      stop("date conversion only implemented for daily, (bi-)weekly and monthly data")
    }
  }
})
setGeneric("epoch<-", function(x, value) standardGeneric("epoch<-"))
setReplaceMethod("epoch", "sts", function(x, value) {
 if (length(value) != nrow(x@observed))
   stop("'epoch' must be of length 'nrow(observed)'")
 if (inherits(value, "Date")) {
   value <- as.integer(value)
   x@epochAsDate <- TRUE
 }
 x@epoch <- value
 x
})
# observed slot
setGeneric("observed", function(x) standardGeneric("observed"))
setMethod("observed", "sts", function(x) {
  return(x@observed)
})
setGeneric("observed<-", function(x, value) standardGeneric("observed<-"))
setReplaceMethod("observed", "sts", function(x, value) {
 x@observed <- value
 x
})
# alarms slot
setGeneric("alarms", function(x) standardGeneric("alarms"))
setMethod("alarms", "sts", function(x) {
  return(x@alarm)
})
setGeneric("alarms<-", function(x, value) standardGeneric("alarms<-"))
setReplaceMethod("alarms", "sts", function(x, value) {
 x@alarm <- value
 x
})
# upperbound slot
setGeneric("upperbound", function(x) standardGeneric("upperbound"))
setMethod("upperbound", "sts", function(x) {
  return(x@upperbound)
})
setGeneric("upperbound<-", function(x, value) standardGeneric("upperbound<-"))
setReplaceMethod("upperbound", "sts", function(x, value) {
 x@upperbound <- value
 x
})
# population slot (actually its populationFrac)
setGeneric("population", function(x) standardGeneric("population"))
setMethod("population", "sts", function(x) {
  return(x@populationFrac)
})
setGeneric("population<-", function(x, value) standardGeneric("population<-"))
setReplaceMethod("population", "sts", function(x, value) {
 x@populationFrac <- value
 x
})
##control slot
setGeneric("control", function(x) standardGeneric("control"))
setMethod("control", "sts", function(x) {
  return(x@control)
})
setGeneric("control<-", function(x, value) standardGeneric("control<-"))
setReplaceMethod("control", "sts", function(x, value) {
 x@control <- value
 x
})
###multinomial Time series slot
##control slot
setGeneric("multinomialTS", function(x) standardGeneric("multinomialTS"))
setMethod("multinomialTS", "sts", function(x) {
  return(x@multinomialTS)
})
setGeneric("multinomialTS<-", function(x, value) standardGeneric("multinomialTS<-"))
setReplaceMethod("multinomialTS", "sts", function(x, value) {
 x@multinomialTS <- value
 x
})

### neighbourhood matrix slot
setGeneric("neighbourhood", function(x) standardGeneric("neighbourhood"))
setMethod("neighbourhood", "sts", function(x) {
  return(x@neighbourhood)
})
setGeneric("neighbourhood<-", function(x, value) standardGeneric("neighbourhood<-"))
setReplaceMethod("neighbourhood", "sts", function(x, value) {
 x@neighbourhood <- value
 x
})


######################################################################
# Miscellaneous access methods
######################################################################

setGeneric("epochInYear", function(x, ...) standardGeneric("epochInYear"))
setGeneric("year", function(x, ...) standardGeneric("year"))


######################################################################
# For stsNC class
######################################################################

### access function for repotringTriangle slot
setGeneric("reportingTriangle", function(x) standardGeneric("reportingTriangle"))
setMethod("reportingTriangle", "stsNC", function(x) {
  return(x@reportingTriangle)
})

### access function for delayCDF slot
setGeneric("delayCDF", function(x) standardGeneric("delayCDF"))
setMethod("delayCDF", "stsNC", function(x) {
  return(x@delayCDF)
})

### access function for SR slot
setGeneric("score", function(x) standardGeneric("score"))
setMethod("score", "stsNC", function(x) {
  return(x@SR)
})

### access function for prediction interval slot
setGeneric("predint", function(x) standardGeneric("predint"))
setMethod("predint", "stsNC", function(x) {
  return(x@pi)
})
