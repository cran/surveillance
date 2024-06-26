################################################################################
### Data structure for CONTINUOUS SPATIO-temporal infectious disease case data
### and a spatio-temporal grid of endemic covariates
###
### Copyright (C) 2009-2018,2021,2024 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
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
#   The column names ".obsInfLength", ".bdist", ".influenceRegion", and ".sources" are reserved.
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
# clipper: engine to use for computing polygon intersections.
######################################################################

obligColsNames_events <- c("time", "tile", "type", "eps.t", "eps.s")
obligColsNames_stgrid <- c("start", "stop", "tile", "area")
dotNames_events <- c(".obsInfLength", ".sources", ".bdist", ".influenceRegion")
reservedColsNames_events <- c(dotNames_events, "BLOCK", "start")
reservedColsNames_stgrid <- c("BLOCK")

as.epidataCS <- function (events, stgrid, W, qmatrix = diag(nTypes),
                          nCircle2Poly = 32, T = NULL,
                          clipper = "polyclip",
                          verbose = interactive())
{
    clipper <- match.arg(clipper)

    # Check and SORT events
    if (verbose) cat("\nChecking 'events':\n")
    events <- check_events(events, verbose = verbose)

    # Check and SORT stgrid
    if (verbose) cat("Checking 'stgrid':\n")
    tiles <- NULL                       # FIXME: add argument to as.epidataCS
    stgrid <- if (missing(stgrid) && inherits(tiles, "SpatialPolygons")) {
        if (verbose) cat("\t(missing, using time-constant 'tiles' grid)\n")
        check_stgrid(tiles2stgrid(tiles, start=0, T=T), verbose = FALSE)
    } else {
        check_stgrid(stgrid, T, verbose = verbose)
    }
    T <- tail(stgrid$stop, 1L)

    # Check class of W and consistency of area
    if (verbose) cat("Checking 'W' ...\n")
    W <- check_W(W, area.other =
                 sum(stgrid[["area"]][seq_len(nlevels(stgrid$tile))]),
                 other = "stgrid")
    stopifnot(identicalCRS(W, events))

    # Check qmatrix
    if (verbose) cat("Checking 'qmatrix' ...\n")
    typeNames <- levels(events$type)
    nTypes <- length(typeNames)     # default value of qmatrix depends on nTypes
    qmatrix <- checkQ(qmatrix, typeNames)

    # Check nCircle2Poly
    stopifnot(isScalar(nCircle2Poly))
    nCircle2Poly <- as.integer(nCircle2Poly)

    # Small helper function converting event index to (time, tile, type) string
    eventidx2string <- function (eventIdx) {
        with(events@data,
             paste(c("time", "tile", "type"), "=",
                   c(time[eventIdx], dQuote(tile[eventIdx]),
                     dQuote(type[eventIdx])),
                   collapse = ", "))
    }

    # Check that all events are part of W
    if (verbose) cat("Checking if all events are part of 'W' ...\n")
    WIdxOfEvents <- over(events, W)
    if (eventNotInWidx <- match(NA, WIdxOfEvents, nomatch = 0L)) {
        stop("the event at (", eventidx2string(eventNotInWidx), ") is not ",
             "inside 'W'")
    }

    # Attach spatio-temporal grid data to events
    events@data <- merge_stgrid(events@data, stgrid, verbose = verbose)
    
    # Calculate observed infection length for log-likelihood
    events$.obsInfLength <- pmin(T - events$time, events$eps.t)

    # Determine potential source events (infective individuals) of each event
    if (verbose) cat("Determining potential event sources ...\n")
    events$.sources <- determineSources(
        eventTimes = events$time, eps.t = events$eps.t,
        eventCoords = coordinates(events), eps.s = events$eps.s,
        eventTypes = events$type, qmatrix = qmatrix)

    # Calculate minimal distance of event locations from the polygonal boundary
    if (verbose) cat("Calculating the events' distances to the boundary ...\n")
    Wowin <- SpP2owin(W)
    events$.bdist <- bdist(coordinates(events), Wowin)

    # Construct spatial influence regions around events
    if (verbose) cat("Constructing spatial influence regions around events ...\n")
    events$.influenceRegion <-
        .influenceRegions(events, Wowin, nCircle2Poly, clipper=clipper)

    # Return components in a list of class "epidataCS"
    res <- list(events = events, stgrid = stgrid, W = W, qmatrix = qmatrix)
    class(res) <- "epidataCS"
    if (verbose) cat("Done.\n\n")
    return(res)
}






######################################################################
# HELPER FUNCTIONS FOR as.epidataCS
######################################################################


### CHECK FUNCTION FOR events ARGUMENT IN as.epidataCS

check_events <- function (events, dropTypes = TRUE, verbose = TRUE)
{
    # Check class and spatial dimensions
    stopifnot(inherits(events, "SpatialPointsDataFrame"))
    if (ncol(events@coords) != 2L) {
        stop("only two spatial dimensions are supported")
    }

    # check suitability of Euclidean geometry
    if (identical(FALSE, is.projected(events))) { # is.projected may return NA
        warning("\"epidataCS\" expects planar coordinates; see 'spTransform'")
    }

    # Check existence of type column
    if (verbose) cat("\tChecking 'type' column ... ")
    events$type <- if ("type" %in% names(events)) {
                       if (dropTypes) factor(events$type) else as.factor(events$type)
                   } else {
                       if (verbose) cat("Setting 'type' to 1 for all events.")
                       factor(rep.int(1L,nrow(events@coords)))
                     }
    if (verbose) cat("\n")

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
    if (verbose) cat("\tChecking 'eps.t' and 'eps.s' columns ...\n")
    with(events@data, stopifnot(is.numeric(eps.t), eps.t > 0,
                                is.numeric(eps.s), eps.s > 0))

    # Transform time into a numeric variable
    if (verbose) cat("\tConverting event time into a numeric variable ...\n")
    events$time <- as.numeric(events$time)
    stopifnot(!is.na(events$time))

    # Check event times for ties
    if (verbose) cat("\tChecking event times for ties ...\n")
    timeIsDuplicated <- duplicated(events$time)
    if (any(timeIsDuplicated)) {
        duplicatedTimes <- sort.int(unique(events$time[timeIsDuplicated]))
        warning("detected concurrent events at ", length(duplicatedTimes),
                " time point", if (length(duplicatedTimes) > 1L) "s", ": ",
                paste(head(duplicatedTimes, 6L), collapse = ", "),
                if (length(duplicatedTimes) > 6L) ", ...")
    }

    # Sort events chronologically
    if (verbose) cat("\tSorting events ...\n")
    events <- events[order(events$time),]

    # First obligatory columns then remainders (epidemic covariates)
    obligColsIdx <- match(obligColsNames_events, names(events@data))
    covarColsIdx <- setdiff(seq_along(events@data), obligColsIdx)
    events@data <- events@data[c(obligColsIdx, covarColsIdx)]
    events@coords.nrs <- numeric(0L)  # forget index of coordinate columns

    # Done.
    return(events)
}



### CHECK FUNCTION FOR stgrid ARGUMENT IN as.epidataCS

check_stgrid <- function (stgrid, T, verbose = TRUE, warn = TRUE)
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
    if (warn && length(reservedCols <- intersect(reservedColsNames_stgrid, names(stgrid)))) {
        warning("replacing existing columns in 'stgrid' which have reserved names: ",
                paste0("'", reservedCols, "'", collapse=", "))
    }

    # Transform tile into a factor variable
    # (also removing unused levels if it was a factor)
    if (verbose) cat("\tConverting 'tile' into a factor variable ...\n")
    stgrid$tile <- factor(stgrid$tile)

    # Transform start times and area into numeric variables
    stgrid$start <- as.numeric(stgrid$start)
    stgrid$area <- as.numeric(stgrid$area)

    # Check stop times
    stgrid$stop <- if (autostop) {
        # auto-generate stop times from start times and T
        if (verbose) cat("\tAuto-generating 'stop' column ...\n")
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
        if (verbose) cat("\tChecking start/stop consistency ...\n")
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
    if (verbose) cat("\tChecking if the grid is complete ...\n")
    blocksizes <- table(stgrid$BLOCK)
    tiletable <- table(stgrid$tile)
    if (length(unique(blocksizes)) > 1L || length(unique(tiletable)) > 1L) {
        stop("'stgrid' is not a full grid")
    }

    # First column BLOCK, then obligCols, then remainders (endemic covariates)
    if (verbose) cat("\tSorting the grid by time and tile ...\n")
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


### MERGE stgrid DATA INTO events@data
## Note: 'events' below refers to the data slot

merge_stgrid <- function (events, stgrid, verbose = TRUE)
{
    # Some basic quantities
    nEvents <- nrow(events)
    timeRange <- with(stgrid, c(start[1], stop[length(stop)]))

    # Are events covered by stgrid?
    if (verbose) {
        cat("Checking if all events are covered by 'stgrid' ...\n")
        ## surveillance > 1.16.0: prehistory events are allowed => BLOCK is NA
        if (events$time[1L] <= timeRange[1L]) {
            cat("  Note: ", sum(events$time <= timeRange[1L]),
                " prehistory events (time <= ", timeRange[1L], ")\n", sep = "")
        }
    }
    if (events$time[nEvents] > timeRange[2L]) {
        stop("found ", sum(events$time > timeRange[2L]),
             " events beyond 'stgrid' (time > ", timeRange[2L], ")")
    }

    # Are all events$tile references really part of the stgrid?
    .events.tile <- factor(events$tile, levels = levels(stgrid$tile))
    if (missingSCellIdx <- match(NA, .events.tile, nomatch = 0L)) {
        stop("the 'events$tile' entry \"", events$tile[missingSCellIdx], "\"",
             " is not a valid level of 'stgrid$tile'")
    }
    events$tile <- .events.tile

    # Map events to corresponding grid cells
    ## FIXME: could use plapply() but then also need a .parallel argument
    if (verbose) cat("Mapping events to 'stgrid' cells ...\n")
    withPB <- verbose && interactive()
    gridcellsOfEvents <- integer(nEvents)
    if (withPB) pb <- txtProgressBar(min=0, max=nEvents, initial=0, style=3)
    for (i in seq_len(nEvents)) {
        gridcellsOfEvents[i] <- gridcellOfEvent(events$time[i], events$tile[i], stgrid)
        if (withPB) setTxtProgressBar(pb, i)
    }
    if (withPB) close(pb)

    # Attach endemic covariates from stgrid to events
    if (verbose) cat("Attaching endemic covariates from 'stgrid' to 'events' ...\n")
    endemicVars <- setdiff(names(stgrid),
                           c(reservedColsNames_stgrid, obligColsNames_stgrid))
    copyCols <- c("BLOCK", "start", endemicVars)
    if (length(replaceCols <- intersect(copyCols, names(events)))) {
        warning("replacing existing columns in 'events' data with ",
                "variables from 'stgrid': ",
                paste0("'", replaceCols, "'", collapse=", "))
        events[replaceCols] <- NULL  # ensure endemic vars are _appended_
    }
    events <- cbind(events, stgrid[gridcellsOfEvents, copyCols])

    return(events)
}



### CHECK FUNCTION FOR W ARGUMENT IN as.epidataCS

check_W <- function (W, area.other = NULL, other, tolerance = 0.001)
{
    W <- as(W, "SpatialPolygons") # i.e. drop data if a SpatialPolygonsDataFrame

    if (!is.null(area.other) && area.other > 0) {
        check_W_area(W, area.other, other, tolerance)
    }

    return(W)
}

check_W_area <- function (W, area.other, other, tolerance = 0.001)
{
    area.W <- areaSpatialPolygons(W)
    if (!isTRUE(all.equal.numeric(area.other, area.W, tolerance = tolerance,
                                  check.attributes = FALSE)))
        warning("area of 'W' (", area.W, ") differs from ",
                "total tile area in '", other, "' (", area.other, ")")
}



### CHECK FUNCTION FOR tiles ARGUMENT IN simEpidataCS()

check_tiles <- function (tiles, levels,
                         events = NULL, areas.stgrid = NULL, W = NULL,
                         keep.data = FALSE, tolerance = 0.05)
{
    stopifnot(inherits(tiles, "SpatialPolygons"),
              is.vector(levels, mode="character"))
    tileIDs <- row.names(tiles)

    ## check completeness of tiles
    if (!identical(tileIDs, levels)) {
        if (any(missingtiles <- !levels %in% tileIDs))
            stop(sum(missingtiles), " regions are missing in 'tiles', ",
                 "check 'row.names(tiles)'")
        ## order tiles by levels and drop any extra tiles
        tiles <- tiles[levels, ]
    }

    ## drop data (also for suitable over-method in check_tiles_events)
    .tiles <- as(tiles, "SpatialPolygons")

    ## check tile specification of events and identical projection
    if (!is.null(events)) {
        check_tiles_events(.tiles, events)
    }

    ## check areas
    areas.tiles <- areaSpatialPolygons(tiles, byid = TRUE)
    if (!is.null(areas.stgrid)) {
        check_tiles_areas(areas.tiles, areas.stgrid, tolerance=tolerance)
    }
    if (!is.null(W)) {
        stopifnot(identicalCRS(tiles, W))
        check_W_area(W, area.other=sum(areas.tiles), other="tiles",
                     tolerance=tolerance)
    }

    ## done
    if (keep.data) tiles else .tiles
}

check_tiles_events <- function (tiles, events)
{
    tiles <- as(tiles, "SpatialPolygons") # remove potential data for over()
    stopifnot(inherits(events, "SpatialPointsDataFrame"),
              identicalCRS(tiles, events))
    tileIDs <- row.names(tiles)
    eventIDs <- row.names(events)

    ## get polygon ID's of events (via overlay)
    eventtiles <- tileIDs[over(events, tiles)]

    if (length(which_not_in_tiles <- which(is.na(eventtiles))))
        warning("some of 'events' are not within 'tiles': ",
                paste0("\"", eventIDs[which_not_in_tiles], "\"", collapse=", "))

    if (!is.null(events@data[["tile"]])) {
        which_disagree <- setdiff(
            which(eventtiles != as.character(events$tile)),
            which_not_in_tiles)
        if (length(which_disagree))
            message("'over(events, tiles)' disagrees with 'events$tile' for events ",
                    paste0("\"", eventIDs[which_disagree], "\"", collapse=", "))
    }
    invisible()
}

check_tiles_areas <- function (areas.tiles, areas.stgrid, tolerance = 0.05)
{
    areas_all_equal <- all.equal.numeric(areas.stgrid, areas.tiles,
                                         tolerance = tolerance,
                                         check.attributes = FALSE)
    if (!isTRUE(areas_all_equal))
        warning("tile areas in 'stgrid' differ from areas of 'tiles': ",
                areas_all_equal)
}


### CONSTRUCT SPATIAL INFLUENCE REGIONS AROUND EVENTS

# An influenceRegion is an object of class "owin" with origin
# at the event (over which we have to integrate by a cubature rule)
# An attribute "area" gives the area of the influenceRegion.
# If it is actually a circular influence region, then there is an attribute
# "radius" denoting the radius of the influence region.
# Argument 'W' can be of class "owin" (preferred) or "SpatialPolygons"
.influenceRegions <- function (events, W, npoly, maxExtent = NULL,
                               clipper = "polyclip")
{
    Wowin <- if (inherits(W, "owin")) W else SpP2owin(W)
    if (is.null(maxExtent)) maxExtent <- diameter.owin(Wowin)
    doIntersection <- switch(
        clipper,  # which package to use for polygon intersection
        "polyclip" = function (center, eps)
            intersectPolyCircle.owin(Wowin, center, eps, npoly),
        ## "rgeos" = function (center, eps) SpP2owin(
        ##     intersectPolyCircle.SpatialPolygons(
        ##         as(W, "SpatialPolygons"), center, eps, npoly)),
        stop("unsupported polygon clipping engine: '", clipper, "'")
        )

    eventCoords <- coordinates(events)
    ## FIXME: could use plapply() but then also need a .parallel argument
    res <- mapply(
        function (x, y, eps, bdist) {
            center <- c(x,y)
            ## if eps is very large, the influence region is the whole region of W
            iR <- shift.owin(
                if (eps > maxExtent) Wowin else doIntersection(center, eps),
                -center)
            ## if iR is actually a circle of radius eps, attach eps as attribute
            attr(iR, "area") <- if (eps <= bdist) {
                attr(iR, "radius") <- eps
                pi * eps^2
            } else area.owin(iR)
            iR
        },
        eventCoords[,1], eventCoords[,2], events$eps.s, events$.bdist,
        SIMPLIFY = FALSE, USE.NAMES = FALSE)

    attr(res, "nCircle2Poly") <- npoly
    attr(res, "clipper") <- clipper
    res
}


### CREATE stgrid TEMPLATE FROM tiles

tiles2stgrid <- function (tiles, start, T)
{
    start <- sort.int(unique.default(start))
    stgrid <- expand.grid(tile = row.names(tiles),
                          start = start,
                          KEEP.OUT.ATTRS = FALSE, stringsAsFactors = TRUE)
    cbind(stgrid,
          stop = rep(c(start[-1L], T), each = length(tiles)),
          area = rep(areaSpatialPolygons(tiles, byid = TRUE), length(start)))
}
