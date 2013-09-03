################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Data structure for CONTINUOUS SPATIO-temporal infectious disease case data
### and a spatio-temporal grid of endemic covariates
###
### Copyright (C) 2009-2012 Sebastian Meyer
### $Revision: 616 $
### $Date: 2013-08-08 15:39:52 +0200 (Don, 08 Aug 2013) $
################################################################################


######################################################################
# MAIN GENERATOR FUNCTION FOR epidataCS OBJECTS
# PARAMS:
# events: SpatialPointsDataFrame of cases with obligatory columns
#   time: time point of event
#   tile: reference to spatial unit (tile) in stgrid, where the event is located
#   type: optional type of event (-> marked twinstim). will be converted to a factor variable.
#   eps.t: maximal temporal influence radius (e.g. length of infectious period, time to culling, etc.), may be Inf
#   eps.s: maximal spatial influence radius (e.g. 100 [km]), may be Inf
#   The remaining columns are further marks of the event, e.g. sex, age of infected person (-> epidemic covariates)
#   The column names "ID", ".obsInfLength", ".bdist", ".influenceRegion", and ".sources" are reserved.
#   "ID": unique chronological ID for the events
#   ".obsInfLength": observed length of the infectious period (being part [0,T])
#   ".bdist": minimal distance of the event locations to the boundary
#   ".influenceRegion": object of class "owin", the intersection of W with b(s,eps.s), with origin at s
#   ".sources": potential sources of infection
# stgrid: data.frame with obligatory columns
#   tile: ID of spatial unit (e.g. id of municipality)
#   start, stop: temporal interval
#   area: area of the spatial unit (tile)
#   The remaining columns are endemic covariates.
#   The column name "BLOCK" is reserved (indexing the time intervals of stgrid).
# W: SpatialPolygons. Observation region. Must have same proj4string as events.
# qmatrix: square indicator matrix (0/1 or TRUE/FALSE) for possible transmission between the event types. will be internally converted to logical. Defaults to an independent spread of the event types.
# nCircle2Poly: accuracy (number of edges) of the polygonal approximation of a circle
# T: end of observation period (=last stop time). Must be specified if only the
#    start but not the stop times are supplied in stgrid (-> auto-generation of stop-times).
######################################################################

obligColsNames_events <- c("time", "tile", "type", "eps.t", "eps.s")
obligColsNames_stgrid <- c("start", "stop", "tile", "area")
reservedColsNames_events <- c("ID", ".obsInfLength", ".sources", ".bdist",
                              ".influenceRegion", "BLOCK", "start")
reservedColsNames_stgrid <- c("BLOCK")

as.epidataCS <- function (events, stgrid, W, qmatrix = diag(nTypes),
                          nCircle2Poly = 32, T = NULL)
{
    # Check and SORT events and stgrid
    cat("\nChecking 'events':\n")
    events <- checkEvents(events)
    cat("Checking 'stgrid':\n")
    stgrid <- checkstgrid(stgrid, T)

    # Check class and proj4string of W and consistency of area
    cat("Checking 'W'...\n")
    W <- as(W, "SpatialPolygons")
    stopifnot(identicalCRS(W, events))
    W.area <- sum(sapply(W@polygons, slot, "area"))
    tiles.totalarea <- sum(stgrid$area[stgrid$BLOCK == 1])
    if (abs(W.area - tiles.totalarea) / max(W.area, tiles.totalarea) > 0.005) {
        cat("\tArea of 'W' =", W.area, "\n")
        cat("\tTotal area of the tiles in 'stgrid' =", tiles.totalarea, "\n")
        warning("area of 'W' should be consistent with the ",
                "total area of the tiles in 'stgrid'")
    }
    
    # Check qmatrix
    cat("Checking 'qmatrix'...\n")
    typeNames <- levels(events$type)
    nTypes <- length(typeNames)
    qmatrix <- checkQ(qmatrix, typeNames)

    # Check nCircle2Poly
    stopifnot(isScalar(nCircle2Poly))
    nCircle2Poly <- as.integer(nCircle2Poly)
    
    # Small helper function converting event index to (time, tile) string
    eventidx2string <- function (eventIdx) {
        paste(c("time", "tile", "type"), "=",
              unlist(events@data[eventIdx,c("time","tile","type")]),
              collapse = ", ")
    }

    # Check that all events are part of W
    cat("Checking if all events are part of 'W'...\n")
    WIdxOfEvents <- over(events, W)
    if (eventNotInWidx <- match(NA, WIdxOfEvents, nomatch = 0L)) {
        stop("the event at (", eventidx2string(eventNotInWidx), ") is not ",
             "inside 'W'")
    }

    # Some basic quantities
    eventCoords <- coordinates(events)
    nEvents <- nrow(eventCoords)
    timeRange <- with(stgrid, c(start[1], stop[length(stop)]))

    # Are event times covered by stgrid?
    cat("Checking if all events are covered by 'stgrid'...\n")
    ## FIXME: what about pre-history events? don't need stgrid-data for them
    if (events$time[1] <= timeRange[1] || events$time[nEvents] > timeRange[2]) {
        stop("event times are not covered by 'stgrid': must be in (",
             timeRange[1L],",",timeRange[2L],"]")
    }

    # Are all events$tile references really part of the stgrid?
    .events.tile <- factor(events$tile, levels = levels(stgrid$tile))
    if (missingSCellIdx <- match(NA, .events.tile, nomatch = 0L)) {
        stop("the 'events$tile' entry \"", events$tile[missingSCellIdx], "\"",
             " is not a valid level of 'stgrid$tile'")
    }
    events$tile <- .events.tile

    # Calculate time point of removal, when event is definitely no longer infective
    removalTimes <- events$time + events$eps.t

    # Calculate distance matrix of events
    cat("Calculating euclidean distance matrix of events...\n")
    eventDists <- as.matrix(dist(eventCoords, method = "euclidean"))
    #diag(eventDists) <- Inf   # infinite distance to oneself (no self-infection), not needed

    # Map events to corresponding grid cells
    # Also precalculate possible origins of events (other infected individuals)
    cat("Mapping events to 'stgrid' cells and",
        "determining potential event sources...\n")
    withPB <- interactive()
    gridcellsOfEvents <- integer(nEvents)
    eventSources <- vector(nEvents, mode = "list")
    if (withPB) pb <- txtProgressBar(min=0, max=nEvents, initial=0, style=3)
    for (i in seq_len(nEvents)) {
        idx <- gridcellOfEvent(events$time[i], events$tile[i], stgrid)
        if (is.na(idx)) {
            stop("could not find information for time point ", events$time[i],
                 " and tile \"", events$tile[i], "\" in 'stgrid'")
        }
        gridcellsOfEvents[i] <- idx
        eventSources[[i]] <- determineSources(
            i, events$time, removalTimes, eventDists[i,], events$eps.s, events$type, qmatrix
        )
        if (withPB) setTxtProgressBar(pb, i)
    }
    if (withPB) close(pb)

    # Attach endemic covariates from stgrid to events
    cat("Attaching endemic covariates from 'stgrid' to 'events'...\n")
    stgridIgnoreCols <- match(setdiff(obligColsNames_stgrid, "start"), names(stgrid))
    copyCols <- setdiff(seq_along(stgrid), stgridIgnoreCols)
    reservedColsIdx <- na.omit(match(names(stgrid)[copyCols], names(events@data),
                                     nomatch=NA_integer_))
    if (length(reservedColsIdx) > 0L) {
        warning("in 'events@data', the existing columns with names of endemic ",
                "covariates from 'stgrid' (",
                paste0("'", names(events@data)[reservedColsIdx], "'", collapse=", "),
                ") have been replaced")
        events@data <- events@data[-reservedColsIdx]
    }
    events@data <- cbind(events@data, stgrid[gridcellsOfEvents, copyCols])

    # Calculate observed infection length = min(T-time, eps.t) for use in log-likelihood
    events$.obsInfLength <- with(events@data, pmin(timeRange[2]-time, eps.t))

    # Attach possible eventSources (infective individuals) to events
    events$.sources <- eventSources

    # Calculate minimal distance of event locations from the polygonal boundary
    cat("Calculating (minimal) distances of the events to the boundary...\n")
    Wgpc <- as(W, "gpc.poly")   # coerce-method from rgeos imported via polyCub
    events$.bdist <- bdist(eventCoords, Wgpc)   # this may take a while
    # avoid conversion to gpc.poly by extracting xylist (first vertex repeated)?
    #lapply(unlist(lapply(W@polygons, slot, "Polygons"), recursive=FALSE), coordinates)


    # Construct spatial influence regions around events
    cat("Constructing spatial influence regions around events...\n")
    events$.influenceRegion <- .influenceRegions(events, Wgpc, W,
                                                 npoly=nCircle2Poly)

    # Attach some useful attributes
    res <- list(events = events, stgrid = stgrid, W = W, qmatrix = qmatrix)
    class(res) <- c("epidataCS", "list")

    cat("Done.\n\n")
    return(res)
}






######################################################################
# HELPER FUNCTIONS FOR as.epidataCS
######################################################################


### CHECK FUNCTION FOR events ARGUMENT IN as.epidataCS

checkEvents <- function (events, dropTypes = TRUE)
{
    # Check class and spatial dimensions
    stopifnot(inherits(events, "SpatialPointsDataFrame"))
    if (ncol(events@coords) != 2L) {
        stop("only two spatial dimensions are supported")
    }

    # Check existence of type column
    cat("\tChecking 'type' column... ")
    events$type <- if ("type" %in% names(events)) {
                       if (dropTypes) factor(events$type) else as.factor(events$type)
                   } else {
                       cat("Setting 'type' to 1 for all events.")
                       factor(rep.int(1L,nrow(events@coords)))
                     }
    cat("\n")

    # Check obligatory columns
    obligColsIdx <- match(obligColsNames_events, names(events), nomatch = NA_integer_)
    if (any(obligColsMissing <- is.na(obligColsIdx))) {
        stop("missing obligatory columns in 'events@data': ",
            paste(obligColsNames_events[obligColsMissing], collapse = ", "))
    }

    # Check other columns on reserved names
    reservedColsIdx <- na.omit(match(reservedColsNames_events, names(events),
                                     nomatch=NA_integer_))
    if (length(reservedColsIdx) > 0L) {
        warning("in 'events@data', the existing columns with reserved names (",
                paste0("'", names(events)[reservedColsIdx], "'", collapse=", "),
                ") have been replaced")
        events@data <- events@data[-reservedColsIdx]
    }    

    # Check that influence radii are numeric and positive (also not NA)
    cat("\tChecking 'eps.t' and 'eps.s' columns...\n")
    with(events@data, stopifnot(is.numeric(eps.t), eps.t > 0,
                                is.numeric(eps.s), eps.s > 0))

    # Transform time into a numeric variable
    cat("\tConverting event time into a numeric variable...\n")
    events$time <- as.numeric(events$time)
    stopifnot(!is.na(events$time))

    # Check event times for ties
    cat("\tChecking event times for ties...\n")
    timeIsDuplicated <- duplicated(events$time)
    if (any(timeIsDuplicated)) {
        duplicatedTimes <- unique(events$time[timeIsDuplicated])
        warning("detected non-unique event times: ",
                "concurrent events at time ",
                if (length(duplicatedTimes) == 1L) "point " else "points\n",
                paste(duplicatedTimes, collapse = ", "))
    }

    # Sort events chronologically
    cat("\tSorting events...\n")
    events <- events[order(events$time),]
    
    # Attribute unique IDs to events (running chronologically from 1 to nEvents)
    events$ID <- seq_along(events)

    # Make ID column the first column, then obligatory columns then remainders (epidemic covariates)
    IDcolIdx <- match("ID", names(events))
    obligColsIdx <- match(obligColsNames_events, names(events))
    covarColsIdx <- setdiff(seq_along(events@data), c(IDcolIdx, obligColsIdx))
    events <- events[c(IDcolIdx, obligColsIdx, covarColsIdx)]

    # Done.
    return(events)
}



### CHECK FUNCTION FOR stgrid ARGUMENT IN as.epidataCS

checkstgrid <- function (stgrid, T)
{
    # Check class
    stopifnot(inherits(stgrid, "data.frame"))

    # Check obligatory columns
    autostop <- FALSE
    if (is.null(stgrid[["stop"]])) {
        if (is.null(T)) stop("'T' must be specified for auto-generation ",
                             "of 'stop' column in 'stgrid'")
        stopifnot(isScalar(T))
        autostop <- TRUE
        stgrid$stop <- NA_real_
    }
    obligColsIdx <- match(obligColsNames_stgrid, names(stgrid), nomatch = NA_integer_)
    if (any(obligColsMissing <- is.na(obligColsIdx))) {
        stop("missing obligatory columns in 'stgrid': ",
            paste(obligColsNames_stgrid[obligColsMissing], collapse = ", "))
    }

    # Check other columns on reserved names
    reservedColsIdx <- na.omit(match(reservedColsNames_stgrid, names(stgrid),
                                     nomatch=NA_integer_))
    if (length(reservedColsIdx) > 0L) {
        warning("in 'stgrid', the existing columns with reserved names (",
                paste0("'", names(stgrid)[reservedColsIdx], "'", collapse=", "),
                ") have been replaced")
        stgrid <- stgrid[-reservedColsIdx]
    }

    # Transform tile into a factor variable
    # (also removing unused levels if it was a factor)
    cat("\tConverting 'tile' into a factor variable...\n")
    stgrid$tile <- factor(stgrid$tile)

    # Transform start times and area into numeric variables
    stgrid$start <- as.numeric(stgrid$start)
    stgrid$area <- as.numeric(stgrid$area)        

    # Check stop times
    stgrid$stop <- if (autostop) {
        # auto-generate stop times from start times and T
        cat("\tAuto-generating 'stop' column...\n")
        starts <- sort(unique(stgrid$start))
        if (T <= starts[length(starts)]) {
            stop("'T' must be larger than the last 'start' time in 'stgrid'")
        }
        stops <- c(starts[-1], T)
        stops[match(stgrid$start, starts)]
    } else {
        as.numeric(stgrid$stop)
    }

    # chronological data.frame of unique periods
    histIntervals <- unique(stgrid[c("start", "stop")])
    histIntervals <- histIntervals[order(histIntervals[,1L]),]
    nBlocks <- nrow(histIntervals)

    if (!autostop) {
        # Check start/stop consistency
        cat("\tChecking start/stop consisteny...\n")
        if (any(histIntervals[,2L] <= histIntervals[,1L])) {
            stop("stop times must be greater than start times")
        }
        startStopCheck <- histIntervals[-1L,1L] != histIntervals[-nBlocks,2L]
        if (startStopCheckIdx <- match(TRUE, startStopCheck, nomatch = 0)) {
            stop("inconsistent start/stop times: time intervals not consecutive ",
                 "at stop time ", histIntervals[startStopCheckIdx,2L])
        }
    }

    # Add BLOCK id
    stgrid$BLOCK <- match(stgrid$start, histIntervals[,1L])

    # Check that we have a full BLOCK x tile grid
    cat("\tChecking if the grid is rectangular (all time-space combinations)...\n")
    blocksizes <- table(stgrid$BLOCK)
    tiletable <- table(stgrid$tile)
    if (length(unique(blocksizes)) > 1L || length(unique(tiletable)) > 1L) {
        stop("'stgrid' is not a full grid")
    }

    # Make BLOCK column the first column, then obligatory columns, then remainders (endemic covariates)
    cat("\tSorting the grid by time and tile...\n")
    BLOCKcolIdx <- match("BLOCK", names(stgrid))
    obligColsIdx <- match(obligColsNames_stgrid, names(stgrid))
    covarColsIdx <- setdiff(seq_along(stgrid), c(BLOCKcolIdx, obligColsIdx))
    stgrid <- stgrid[c(BLOCKcolIdx, obligColsIdx, covarColsIdx)]

    # Sort by BLOCK and tile
    stgrid <- stgrid[order(stgrid$BLOCK, stgrid$tile),]

#     # Get row indexes of the blocks' first/last rows
#     beginBlock <- match(seq_len(nBlocks), stgrid[["BLOCK"]])
#     endBlock <- c(beginBlock[-1L]-1L, nrow(stgrid))

    # Done.
    return(stgrid)
}



### CONSTRUCT SPATIAL INFLUENCE REGIONS AROUND EVENTS

# An influenceRegion is an object of class "owin" with origin
# at the event (over which we have to integrate by a cubature rule)
# An attribute "area" gives the area of the influenceRegion.
# If it is actually a circular influence region, then there is an attribute
# "radius" denoting the radius of the influence region.
.influenceRegions <- function (events, Wgpc, W, npoly, maxExtent = NULL)
{
    if (is.null(maxExtent)) maxExtent <- maxExtent.gpc.poly(Wgpc)
    eventCoords <- coordinates(events)
    nEvents <- nrow(eventCoords)
    res <- vector(nEvents, mode = "list")
    for (i in seq_len(nEvents)) {
        eps <- events$eps.s[i]
        center <- eventCoords[i,]
        iRgpc <- if (eps > maxExtent) { # influence region is whole region of W
                Wgpc
            } else { # influence region is a subset of W
                intersectCircle(Wgpc, W, center, eps, npoly)
            }
        res[[i]] <- as(scale(iRgpc, center=center), "owin") # coerce via polyCub
        ## if iR is actually a circle of radius eps, attach eps as attribute
        r <- if (eps <= events$.bdist[i]) eps else NULL
        attr(res[[i]], "radius") <- r
        attr(res[[i]], "area") <- if(is.null(r)) area.owin(res[[i]]) else pi*r^2
    }
    attr(res, "nCircle2Poly") <- npoly
    return(res)
}
