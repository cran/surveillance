################################################################################
### Conversion between "ts" and "sts", and from "sts" to "data.frame"
###
### Copyright (C) 2014 Michael Hoehle, 2015-2017,2019 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################

### Convert a simple "ts" object to an "sts" object

setAs(from = "ts", to = "sts", def = function (from) {
    ## Extract frequency and start from the "ts" object
    freq <- frequency(from)
    start <- start(from)
    if (length(start) == 1)
        stop("could not convert time series start() to (year, index) form")

    ## Remove "tsp" attribute and "ts"/"mts" class
    tsp(from) <- NULL
    ## "tsp<-"(x,NULL) is documented to also remove "ts" and "mts" classes
    ## but in R < 3.3.0, it did not remove "mts" (see PR#16769)
    from <- unclass(from)

    ## Create the sts object
    .sts(observed = from, start = start, freq = freq)
})


### Convert an "sts" object to a simple "ts" object

as.ts.sts <- function (x, ...)
{
    ts(data = x@observed, start = x@start, frequency = x@freq)
}
setAs(from = "sts", to = "ts", def = function (from) as.ts.sts(from))


### Convert an "sts" object to an eXtensible Time Series "xts"

as.xts.sts <- function (x, order.by = epoch(x, as.Date = TRUE), ...)
{
    if (!missing(order.by) || x@freq %in% c(52, 365)) {
        xts::xts(x = x@observed, order.by = order.by, ...)
    } else {
        ## frequencies 4 and 12 are nicely handled by the as.xts.ts method
        xts::as.xts(as.ts.sts(x), ...)
    }
}


### Convert an "sts" object to a data frame suitable for regression

as.data.frame.sts <- function(x, row.names = NULL, optional = FALSE, # from the generic
                              tidy = FALSE, as.Date = x@epochAsDate, ...)
{
  if (tidy)
    return(tidy.sts(x, ...))

  #Convert object to data frame and give names
  res <- data.frame("observed" = x@observed,
                    "epoch" = epoch(x, as.Date = as.Date),
                    "state" = x@state,
                    "alarm" = x@alarm,
                    "upperbound" = x@upperbound,
                    "population" = x@populationFrac,
                    check.names = FALSE)

  names(res) <- if (ncol(x) > 1) {
    ## names from data.frame() above should already be as intended
    namesObs <- colnames(x@observed, do.NULL = FALSE, prefix = "observed")
    c(paste0("observed.", namesObs), "epoch",
      paste0("state.", namesObs), paste0("alarm.", namesObs),
      paste0("upperbound.", namesObs), paste0("population.", namesObs))
  } else {
    c("observed", "epoch", "state", "alarm", "upperbound", "population")
  }

  #Find out how many epochs there are each year
  res$freq <- if (x@epochAsDate) {
    year <- strftime(epoch(x), if (x@freq == 52) "%G" else "%Y")
    epochStr <- switch(as.character(x@freq),
                       "12" = "%m", "52" = "%V", "365" = "%j")
    maxEpoch <- vapply(X = unique(year), FUN = function (Y) {
        dummyDates <- as.Date(paste0(Y, "-12-", 26:31))
        max(as.numeric(strftime(dummyDates, epochStr)))
    }, FUN.VALUE = 0, USE.NAMES = TRUE)
    maxEpoch[year]
  } else { # just replicate the fixed frequency
    x@freq
  }

  #Add a column denoting the epoch fraction within the current year
  res$epochInPeriod <- epochInYear(x) / res$freq

  return(res)
}

setMethod("as.data.frame", signature(x = "sts"), as.data.frame.sts)


### convert an "sts" object to a "data.frame" in long (tidy) format

tidy.sts <- function (x, ...)
{
    unitNames <- colnames(x, do.NULL = FALSE, prefix = "observed")
    v.names <- c("observed", "state", "alarm", "upperbound", "population")
    stswide <- as.data.frame(x, tidy = FALSE, as.Date = FALSE)
    ## nrow(stswide) = nrow(x), i.e., one row per epoch
    stswide$year <- year(x)
    stswide$epochInYear <- epochInYear(x)
    stswide$date <- tryCatch(
        epoch(x, as.Date = TRUE),  # only works for particular values of x@freq
        error = function (e) {message("Note: ", e$message); as.Date(NA)}
    )
    if ((nUnit <- ncol(x)) == 1L) {
        stslong <- data.frame(stswide, "unit" = factor(unitNames),
                              check.names = FALSE)
    } else {
        ## we have observed/population/... columns for each unit
        varying <- sapply(X = v.names, FUN = paste, unitNames, sep = ".",
                          simplify = FALSE, USE.NAMES = TRUE)
        stslong <- reshape(
            data = stswide, direction = "long",
            varying = varying, v.names = v.names,
            timevar = "unit", times = unitNames,
            idvar = "epoch")
        stslong$unit <- factor(stslong$unit, levels = unitNames)
        attr(stslong, "reshapeLong") <- NULL
    }
    row.names(stslong) <- NULL
    ## reorder variables (ordering from above differs depending on nUnit)
    stslong[c("epoch", "unit",
              "year", "freq", "epochInYear", "epochInPeriod", "date",
              v.names)]
}
