################################################################################
### Initialization and other basic methods for the S4 class "sts"
###
### Copyright (C) 2007-2014 Michael Hoehle
### Copyright (C) 2012-2019,2021,2023-2025 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################


######################################################################
# initialize-method -- see ../man/sts-class.Rd for class information
######################################################################

#Ensure that all matrix slots have the same dimnames, which are
#always taken from the observed matrix
fix.dimnames <- function(x) {
  dn <- dimnames(x@observed)
  #Make sure all arrays have the same dimnames
  dimnames(x@alarm) <- dimnames(x@state) <- dimnames(x@upperbound) <-
    dimnames(x@populationFrac) <- dn
  #Special for neighbourhood
  dimnames(x@neighbourhood) <- dn[c(2L,2L)]
  return(x)
}

## a user-level constructor function,
## which calls the standard generator function .sts(),
## which calls initialize() on the "sts" prototype - see init.sts() below
## NOTE: using sts() is the preferred approach since surveillance 1.10-0
## NOTE: NULL arguments are ignored => default slot values
sts <- function (observed,
                 start = c(2000, 1), frequency = 52, # prototype values
                 epoch = NULL, # defaults to 1:nrow(observed), can be Date
                 population = NULL, # an alias for "populationFrac"
                 map = NULL, # allow sf input for the map slot
                 ...) # further named arguments representing "sts" slots
{
    slots <- list(observed = observed, start = start, freq = frequency,
                  epoch = epoch, map = map, ...)

    if (!is.null(population)) {
        if ("populationFrac" %in% names(slots))
            warning("'population' takes precedence over 'populationFrac'")
        slots$populationFrac <- population
    } # else "populationFrac" is a possible element of ...

    if (inherits(epoch, "Date")) {
        ## FIXME: guess missing start value similar to linelist2sts
        ## if (missing(start) && frequency == 52)
        ##     slots$start <- unlist(isoWeekYear(epoch[1L]), use.names = FALSE)
        slots$epoch <- as.integer(epoch)
        slots$epochAsDate <- TRUE
    }

    if (inherits(map, "sf")) {
        if (!requireNamespace("sf", quietly = TRUE)) {
            stop("'map' is of class \"sf\" but package 'sf' is not installed.\n",
                 "  Alternatively, use a map of class \"SpatialPolygons\" (from 'sp').")
        }
        ## allow "sf" input and transform it to "SpatialPolygons*"
        slots$map <- sf::as_Spatial(map)
    } else if (!is.null(map) && !inherits(map, "SpatialPolygons")) {
        stop("'map' must inherit from \"sf\" or \"SpatialPolygons\"")
    }

    ## call the standard generator function with explicitly set slots
    isNULL <- vapply(X = slots, FUN = is.null,
                     FUN.VALUE = FALSE, USE.NAMES = FALSE)
    do.call(.sts, slots[!isNULL])
}

## initialize-method called by new("sts", ...),
## the long-standing default way of creating "sts" objects.
## For backward-compatibility, we keep this customized initialize-method,
## although it would be cleaner to put things into the generator function
## and use the default initialize-method.
init.sts <- function(.Object, ..., # also for slots of classes extending "sts"
                     observed, # use copy constructor if missing(observed)
                     ## the following default arguments depend on dim(observed)
                     epoch = seq_len(nTime),
                     state = matrix(FALSE, nTime, nUnit),
                     alarm = matrix(NA, nTime, nUnit),
                     upperbound = matrix(NA_real_, nTime, nUnit),
                     neighbourhood = matrix(NA, nUnit, nUnit),
                     populationFrac = matrix(1/nUnit, nTime, nUnit),
                     ## FIXME: change default to a matrix of NA_real_ ?
                     ## the map slot needs special treatment (see below)
                     map = .Object@map # old/prototype value
                     ## the remaining slots have useful prototype values
                     ## and are handled as part of ...
                     ##start = c(2000, 1), freq = 52,
                     ##epochAsDate = FALSE, multinomialTS = FALSE,
                     ##control = .Object@control
                     )
{
    if (nargs() < 2) # nothing to do
        return(.Object)

    if (missing(observed)) { # use default initialize-method
        ## such that, e.g., initialize(stsObj, map=newMap) will set a new map
        ## and copy other slots from stsObj instead of (re-)setting to defaults,
        ## as well as to support new("stsBP", stsObj, ci=ci, lambda=lambda).
        ## CAVE: automatic dimension correction of matrix slots is not done.
        .Object <- callNextMethod()
        ## Drawback: .Object@map has been forced to "SpatialPolygons"
        if (!missing(map)) # restore the supplied map
            .Object@map <- map
        ## If missing(map), .Object@map = as(stsObj@map, "SpatialPolygons"),
        ## i.e., data will be lost => map=stsObj@map must be passed explicitly
        .Object <- fix.dimnames(.Object)
        return(.Object)
    }

    ## Ensure matrix form (auto-conversion is useful for single time series)
    observed <- as.matrix(observed)
    nUnit <- ncol(observed)
    nTime <- nrow(observed)
    state <- as.matrix(state)
    alarm <- as.matrix(alarm)
    upperbound <- as.matrix(upperbound)

    ## clear rownames and set colnames for the matrix of observed counts
    if (is.null(namesObs <- colnames(observed))){
        namesObs <- paste0("observed", seq_len(nUnit))
    }
    dimnames(observed) <- list(NULL, namesObs)

    ## if there is only one state-vector for more than one area, repeat it
    if (nUnit > 1 && ncol(state) == 1 && length(state) == nTime) {
        state <- rep.int(state, nUnit)
        dim(state) <- c(nTime, nUnit)
    }

    ## time-constant population fractions can be provided as a single vector
    if (is.vector(populationFrac, mode = "numeric")) {
        if (length(populationFrac) != nUnit)
            stop("population vector has length ", length(populationFrac),
                 " but ", nUnit, " units are 'observed'", call. = FALSE)
        populationFrac <- matrix(populationFrac, nTime, nUnit, byrow=TRUE)
    }

    ## we need to set the map manually since the initialize,ANY-method called
    ## next would coerce a "SpatialPolygonsDataFrame" to "SpatialPolygons"
    if (!missing(map))
        .Object@map <- map

    ## set all other slots (including for classes extending this class)
    ## using the default initialize-method
    .Object <- callNextMethod(.Object, ...,
                              observed=observed, epoch=epoch,
                              state=state, alarm=alarm, upperbound=upperbound,
                              neighbourhood=neighbourhood,
                              populationFrac=populationFrac)
    ## this also checks validObject(.Object)
    ## for nUnit > 1, it will catch if any column names differ from namesObs

    ## use dimnames(observed) for all matrix slots (could be unnamed)
    .Object <- fix.dimnames(.Object)

    return(.Object)
}

setMethod("initialize", "sts", init.sts)


###########################################################################
# Conversion between old "disProg" and new "sts" classes
###########################################################################

## transform a "disProg" object to the new "sts" class
disProg2sts <- function(disProgObj, map=NULL) {
  disProgObj$map <- map
  ## NOTE: we cannot trust disProgObj$week to be a valid "epoch" specification,
  ## e.g., the week in data("ha") refers to the week number _within_ a year.
  ## CAVE: in "disProg" objects, several elements may be missing or NULL,
  ## and there could be further elements not matching any "sts" slot,
  ## e.g., in "disProg" objects generated by sim.pointSource()
  validElements <- names(disProgObj) %in% slotNames("sts") &
      !vapply(X=disProgObj, FUN=is.null, FUN.VALUE=FALSE, USE.NAMES=FALSE)
  ## initialize an "sts" object using the valid "disProg" elements
  stsObj <- do.call(.sts, disProgObj[validElements])
  return(stsObj)
}

## The reverse action
sts2disProg <- function(sts) {
  disProgObj <- list(  # no need for create.disProg() checks, sts is formal
    week=sts@epoch, observed=sts@observed, state=sts@state,
    start=sts@start, freq=sts@freq,
    neighbourhood=sts@neighbourhood, populationFrac=sts@populationFrac,
    epochAsDate=sts@epochAsDate
  )
  class(disProgObj) <- "disProg"
  #For survRes: alarm=sts@alarm, upperbound=sts@upperbound)
  return(disProgObj)
}


###########################################################################
#Method to aggregate over all units, either the time series is aggregated
#so a new sampling frequency of nfreq units per time slot is obtained.
#The other alternative is to aggregate all units.
#
# Note: The function is not 100% consistent with what the generic
#       aggregate does.
#
# Warning: In case the aggregation is by unit the upperbound slot is set
#          to NA. Furthermore the MAP object is left as.is, but
#          the object cannot be plotted anymore.
#
# Params:
#   by - a string being either "time" or "unit"
#   nfreq - new sampling frequency if by=="time". If "all" then all
#           time instances are summed.
###########################################################################

aggregate.sts <- function(x, by="time", nfreq="all", ...)
{
  by <- match.arg(by, choices = c("time", "unit"))

  #Aggregate time
  if (by == "time") {
    if (nfreq == "all") {
      howmany <- dim(x@observed)[1]
    } else if (nfreq == x@freq) { # nothing to do
      return(x)
    } else { # nfreq != x@freq
      howmany <- x@freq / nfreq
      if (howmany - ceiling(howmany) != 0)
          stop("nfreq has to be a multiple of x@freq.")
    }

    n <- dim(x@observed)[1]
    m <- ceiling(n/howmany)
    new <- rep(1:m,each=howmany)[1:n]
    x@freq <- ifelse(nfreq == "all", howmany, nfreq)
    x@epoch <- 1:m

    x@observed <- as.matrix(aggregate(x@observed,by=list(new),sum)[,-1])
    x@state <- as.matrix(aggregate(x@state,by=list(new),sum)[,-1])>0
    x@alarm <- as.matrix(aggregate(x@alarm,by=list(new),sum)[,-1]) # number of alarms
    x@upperbound <- as.matrix(aggregate(x@upperbound,by=list(new),sum)[,-1])

    ## summing population (fractions) over time
    had_fractions <- !x@multinomialTS && all(rowSums(x@populationFrac) == 1)
    x@populationFrac <- as.matrix(aggregate(x@populationFrac,by=list(new),sum)[,-1])
    if (isTRUE(had_fractions)) { # population fractions need to be recomputed
      x@populationFrac <- x@populationFrac / rowSums(x@populationFrac)
    }
  }

  #Aggregate units
  if (by == "unit") {
    #Aggregate units
    x@observed <- as.matrix(rowSums(x@observed))
    x@state <- as.matrix(rowSums(x@state))>0
    x@alarm <- as.matrix(rowSums(x@alarm))>0 # contrary to counting for by="time"!
    #There is no clever way to aggregate the upperbounds
    x@upperbound <- matrix(NA_real_,ncol=ncol(x@alarm),nrow=nrow(x@alarm))
    x@populationFrac <- as.matrix(rowSums(x@populationFrac))
    x@neighbourhood <- matrix(NA, 1, 1) # consistent with default for new("sts")
    ## we have lost colnames
    colnames(x@observed) <- "overall"
    x <- fix.dimnames(x)
    ## drop the map (set to empty prototype)
    x@map <- new(getSlots("sts")[["map"]])
  }

  #validObject(x) #just a check

  return(x)
}

setMethod("aggregate", signature(x="sts"), aggregate.sts)


#####################################################################
# Miscellaneous access methods
####################################################################

setMethod("dim", "sts", function (x) dim(x@observed))
setMethod("dimnames", "sts", function (x) dimnames(x@observed))
setMethod("frequency", "sts", function (x, ...) x@freq)
setMethod("start", "sts", function (x, ...) x@start)

## ## time() method to extract fractional year index (not really needed)
## setMethod("time", "sts", function (x, ...) {
##     tsp1 <- x@start[1L] + (x@start[2L] - 1)/x@freq
##     seq.int(tsp1, by = 1/x@freq, length.out = length(x@epoch))
## })

## Extract which observation within year we have
## (this could have been named "cycle", for consistency with "ts")
setMethod("epochInYear", "sts", function(x,...) {
  if (x@epochAsDate && x@freq %in% c(12, 52, 365)) {
    epochStr <- switch(as.character(x@freq),
                       "12" = "%m", "52" = "%V", "365" = "%j")
    as.numeric(strftime(epoch(x), epochStr))
  } else {
    index <- if (x@epochAsDate) { # non-standard frequency
                 seq_along(x@epoch)
             } else x@epoch  # should always be 1:nrow(x) actually
    (index-1 + x@start[2]-1) %% x@freq + 1
  }
})

#Extract the corresponding year for each observation
setMethod("year", "sts", function(x,...) {
  if (x@epochAsDate) {
    as.numeric(strftime(epoch(x), if (x@freq == 52) "%G" else "%Y"))
  } else {
    ((x@epoch-1 + x@start[2]-1) + (x@freq*x@start[1])) %/% x@freq
  }
})


#####################################################################
#[-method for truncating the time series and/or selecting units
#####################################################################

setMethod("[", "sts", function(x, i, j, ..., drop = FALSE) {
  nTimeOriginal <- nrow(x@observed)
  if (missing(i)) { # set default value
    i <- seq_len(nTimeOriginal)
  } else if (anyNA(i)) {
    stop("missing row index values are not supported")
  } else if (is.logical(i)) { # convert to integer index
    i <- which(rep_len(i, nTimeOriginal))
  } else if (is.character(i)) {
    stop("character row indices are not supported")
  } else if (any(i < 0)) { # convert to (positive) indices
    if (any(i > 0)) stop("only 0's may be mixed with negative subscripts")
    i <- setdiff(seq_len(nTimeOriginal), -i)
  } else if (any(i0 <- i == 0)) { # drop 0's (for the diff check below)
    i <- i[!i0]
  }
  ## if(missing(j)) j <- seq_len(ncol(x@observed))   # redundant
  if (!missing(j) && anyNA(j))
    stop("missing column index values are not supported")
  ## FIXME: should probably warn about duplicated column indices

  ## check if i is a regular integer sequence (not invalidating freq)
  if (any(diff(i) != 1))
    warning("irregular row index could invalidate \"freq\"")

  x@epoch <- x@epoch[i]
  x@observed <- x@observed[i,j,drop=FALSE]
  x@state <- x@state[i,j,drop=FALSE]
  x@alarm <- x@alarm[i,j,drop=FALSE]

  recompute_fractions <- !missing(j) && !x@multinomialTS &&
      all(rowSums(x@populationFrac) == 1)
  x@populationFrac <- x@populationFrac[i,j,drop=FALSE]
  if (isTRUE(recompute_fractions)) {
    x@populationFrac <- x@populationFrac / rowSums(x@populationFrac)
   }
  x@upperbound <- x@upperbound[i,j,drop=FALSE]

  #Neighbourhood matrix
  if (ncol(x@observed) != ncol(x@neighbourhood) &&  # selected units
      !all(x@neighbourhood %in% c(NA,0,1))) { # no adjacency matrix
      message("Note: selection of units could invalidate the 'neighbourhood'")
      ## e.g., if 'neighbourhood' specifies neighbourhood orders
  }
  x@neighbourhood <- x@neighbourhood[j,j,drop=FALSE]

  #Fix the "start" and "epoch" entries (if necessary)
  if (any(i != 0) && i[1] != 1) {
    #Note: This code does not work if we have week 53s!
    i.min <- min(i)  # in regular use, this should actually be i[1]
    new.sampleNo <- x@start[2] + i.min - 1
    start.year <- x@start[1] + (new.sampleNo - 1) %/% x@freq
    start.sampleNo <- (new.sampleNo - 1) %% x@freq + 1
    x@start <- c(start.year, start.sampleNo)
    if (!x@epochAsDate) {
      ## we also have to update epoch since it is relative to start
      ## and actually it should always equal 1:nrow(observed)
      x@epoch <- x@epoch - i.min + 1L
    }
    ## if (x@epochAsDate && x@freq == 52) {
    ##   ## FIXME: should we derive start from the first date?
    ##   ISO <- isoWeekYear(as.Date(x@epoch[1], origin = "1970-01-01"))
    ##   x@start <- c(ISO$ISOYear, ISO$ISOWeek)
    ## }
  }

  ## Note: We do not automatically subset the map according to j, since
  ##       identical(row.names(map), colnames(observed))
  ##       is not a property of the sts-class; Unmonitored regions are allowed.
  ##       The map can also be empty (prototype value).
  if (drop && !missing(j)) {
      if (!is.character(j))
          stop("'drop = TRUE' requires character-type column indices")
      if (length(x@map))
          x@map <- x@map[j,]
      else
          warning("nothing to drop; object has no map")
  }

  #Done
  return(x)
})


#########################################################################
## Plot method ... the type argument specifies what type of plot to make
##
## plot as multivariate time series:  type = observed ~ time | unit
## plot overall time series:          type = observed ~ time
## plot as map aggregated over time:  type = observed ~ unit
## the specific plot functions are in separate files (stsplot_*.R)
########################################################################

plot.sts <- function (x, type = observed ~ time | unit, ...)
{
    ## Valid formula?
    stopifnot(inherits(type, "formula"))
    LHS <- if ((length(type) == 2) || type[[2]] == "observed") "observed"
           else if (type[[2]] == "alarm") "alarm"
           else stop("'type' must have LHS 'observed', 'alarm' or empty")
    RHS <- type[[length(type)]]
    RHS_elts <- all.names(RHS)
    valid <- RHS_elts %in% c("1","unit","|","time")  # 1 is deprecated
    if (!all(valid))
        stop("not a valid plot 'type'")

    if ("time" %in% RHS_elts) {
        ## Temporal plots
        if (LHS == "alarm") {
            ## no auto-aggregation: alarm~time == alarm~time|unit
            stsplot_alarm(x, ...)
        } else {
            if (RHS == "time" && ncol(x) > 1)
                x <- aggregate(x, by = "unit")
            stsplot_time(x, ...)
        }
    } else {
        ## Spatial plots
        if (RHS != "unit") # was stsplot_spacetime(x, type, ...)
            warning(sprintf("plot type '%s' is defunct; using '%s'",
                            paste(deparse(type), collapse = ""),
                            "observed ~ unit"), call. = FALSE)
        stsplot_space(x, ...)
    }
}

setMethod("plot", signature(x="sts", y="missing"), plot.sts)


## define how "sts" objects get printed

setMethod( "show", "sts", function( object ){
  cat( "-- An object of class ", class(object), " -- \n", sep = "" )
  if (!object@epochAsDate) {
    cat( "freq:\t\t", object@freq,"\n" )
  } else {
    epochStr <- switch( as.character(object@freq), "12" = "%m","52" =  "%V","365" = "%j")
    cat( "freq:\t\t", paste(object@freq," with strptime format string ",epochStr,"\n",sep=""))
  }
  if (!object@epochAsDate) {
    cat( "start:\t\t",object@start,"\n" )
  } else {
    cat( "start:\t\t",paste(epoch(object)[1]),"\n" )
  }
  cat( "dim(observed):\t", dim(object@observed), "\n\n")

  n <- 1
  cat("Head of observed:\n")
  print(head(object@observed,n))

  if (npoly <- length(object@map)) {
      cat("\nmap:", npoly, "Polygons, ")
      ## no longer print SpatialPolygons as this may load heavy sf (for is.projected)
      ## print(modifyList(summary(object@map), list(data=NULL))) # no data summary
      if (inherits(object@map, "SpatialPolygonsDataFrame")) {
          cat(ncol(object@map), "variables")
          cat("\nHead of map@data:\n")
          print(head(object@map@data, n))
      } else cat("without data\n")
  }

  if (ncol(object@observed) > 1 && !all(is.na(object@neighbourhood))) {
      cat("\nHead of neighbourhood:\n")
      print( head(object@neighbourhood,n))
  }
} )
