################################################################################
# Author: Sebastian Meyer
# Date: 19 Jun 2009
#
# This file contains methods related to the class "epidata" (a history
# data.frame of the observed epidemic): A converter function as.epidata and 
# standard methods like print, summary and plot.
# There are also more sophisticated functions 'animate' and 'stateplot'.
# NB: The summary and plot functions stem from my Bachelor Thesis.
################################################################################


################################################################################
# CONVERTER FUNCTION
# the main purpose is to prepare the data of an epidemic for the inference
# function 'twinSIR', i.e. perform consistency checks, transform data to common
# format and calculate epidemic covariates with distance function f.
# - we assume fixed coordinates (this is important as time-varying coordinates
# would result in more sophisticated and time consuming calculations of
# distance matrices) !
# - in the first block (start = t0) all id's must be present (for coordinates)
# - those id's with atRiskY(t0) = 0 are taken as initially infectious
# - SIS epidemics are possible, but must be given as SIRS with pseudo R-events,
# i.e. individuals will be removed and become susceptible directly afterwards
################################################################################

# default method (original data in a data.frame)
as.epidata.default <- function(data, id.col, start.col, stop.col, atRiskY.col,
    event.col, Revent.col, coords.cols, f = list(), ...)
{
    cl <- match.call()
    
    # If necessary, convert 'data' into a data.frame (also converting
    # column names to syntactically correct names for use in formulae)
    data <- as.data.frame(data, stringsAsFactors = FALSE)
    
    # Check f
    if (length(f) > 0L) {
        if (is.null(coords.cols)) {
            stop("need coordinates to calculate epidemic covariates")
        }
        if (!is.list(f) || is.null(names(f)) || any(!sapply(f, is.function))) {
            stop("'f' must be a named list of functions or 'list()'")
        }
        if (any(fInfNot0 <- sapply(f, function(B) B(Inf)) != 0)) {
            stop("all functions in 'f' must return 0 at infinite distance: ",
                 "f[[i]](Inf) == 0 fails for i = ",
                 paste(which(fInfNot0), collapse=", "))
        }
    }
    
    # Use column numbers as indices and check them
    colargs <- c("id.col", "start.col", "stop.col", "atRiskY.col",
                 "event.col", "Revent.col", "coords.cols")
    colidxs <- structure(as.list(numeric(length(colargs))), names = colargs)
    for (colarg in colargs) {
        colidx <- get(colarg, inherits = FALSE)
        if (colarg != "coords.cols" && length(colidx) != 1L) {
            stop("the column specifier '", colarg, "' must be of length 1")
        }
        if (is.character(colidx)) {
            colidx <- match(colidx, colnames(data))
            if (any(is.na(colidx))) {
                stop("'", colarg, " = ", deparse(cl[[colarg]]), "': ",
                     "column does not exist in 'data'")
            }
        } else if (is.numeric(colidx) && any(colidx<1L | colidx>ncol(data))) {
            stop("'", colarg, " = ", deparse(cl[[colarg]]), "': ",
                 "column index must be in [1; ", ncol(data), "=ncol(data)]")
        }
        colidxs[[colarg]] <- colidx
    }
    
    # Rename main columns to default column names
    colidxsVec <- unlist(colidxs)
    colnams <- c("id", "start", "stop", "atRiskY", "event", "Revent")
    colnames(data)[colidxsVec[1:6]] <- colnams
    usedReservedName <- any(colnams %in% colnames(data)[-colidxsVec[1:6]])
    
    # REORDER COLUMNS, so that main columns come first (also for make.unique)
    data <- data[c(colidxsVec, setdiff(seq_len(NCOL(data)), colidxsVec))]
    
    # Make columns names unique (necessary if other column with name in colnams)
    if (usedReservedName) {
        colnames(data) <- make.unique(colnames(data))
        message("Some other columns had reserved names and have been renamed")
    }
    
    # Convert id into a factor (also removing unused levels if it was a factor)
    data[["id"]] <- factor(data[["id"]])
    
    # Check atRiskY, event and Revent for values other than 0 and 1
    for (var in c("atRiskY", "event", "Revent")) {
        data[[var]] <- as.numeric(data[[var]])
        if (any(! data[[var]] %in% c(0,1)))
            stop("'", var, "' column may only assume values 0 and 1")
    }
    
    # Check consistency of atRiskY and event (event only if at-risk)
    noRiskButEvent <- data[["atRiskY"]] == 0 & data[["event"]] == 1
    if (noRiskButEventRow <- match(TRUE, noRiskButEvent, nomatch = 0)) {
        stop("inconsistent atRiskY/event indicators in row ",
             noRiskButEventRow, ": event only if at risk")
    }

    # Check event (infection) times for ties
    eventTimes <- data[data[["event"]] == 1, "stop"]
    ReventTimes <- data[data[["Revent"]] == 1, "stop"]
    duplicatedEventTime <- duplicated(c(eventTimes, ReventTimes))
    if (duplicatedEventTimeIdx <- match(TRUE, duplicatedEventTime, nomatch=0)) {
        stop("non-unique event times: concurrent event/Revent at time ",
             c(eventTimes, ReventTimes)[duplicatedEventTimeIdx])
    }

    # Check start/stop consistency and add block id
    histIntervals <- unique(data[c("start", "stop")])
    histIntervals <- histIntervals[order(histIntervals[,1L]),]
    nBlocks <- nrow(histIntervals)
    if (any(histIntervals[,2L] <= histIntervals[,1L])) {
        stop("stop times must be greater than start times")
    }
    startStopCheck <- histIntervals[-1L,1L] != histIntervals[-nBlocks,2L]
    if (startStopCheckIdx <- match(TRUE, startStopCheck, nomatch = 0)) {
        stop("inconsistent start/stop times: time intervals not consecutive ",
             "at stop time ", histIntervals[startStopCheckIdx,2L])
    }
    if ("BLOCK" %in% colnames(data)) {
        warning("column name 'BLOCK' is reserved, ",
                "existing column has been replaced")
    }
    data[["BLOCK"]] <- match(data[["start"]], histIntervals[,1L])
    
    # SORT by block/id and create indexes for block borders
    data <- data[order(data[["BLOCK"]], data[["id"]]),]
    beginBlock <- match(seq_len(nBlocks), data[["BLOCK"]])
    endBlock <- c(beginBlock[-1L]-1L, nrow(data))
    
    # make block column the first column
    BLOCK.col <- match("BLOCK", colnames(data))
    data <- data[c(BLOCK.col, setdiff(seq_along(data), BLOCK.col))]
    coords.cols <- 1L + 6L + seq_along(colidxs[["coords.cols"]])

    # Check consistency of atRiskY and event (not at-risk after event) 
    .checkFunction <- function(eventblock, eventid)
    {
        rowsOfNextBlock <- beginBlock[eventblock+1L]:endBlock[eventblock+1L]
        nextBlockData <- data[rowsOfNextBlock, c("id", "atRiskY")]
        idIdx <- which(nextBlockData[["id"]] == eventid)
        if (length(idIdx) == 1L && nextBlockData[idIdx, "atRiskY"] == 1) {
            stop("inconsistent atRiskY/event indicators for id '", eventid,
                 "': should not be at risk immediately after event")
        }
    }
    eventTable <- data[data[["event"]] == 1,]
    for(k in seq_len(nrow(eventTable)))
    {
        .checkFunction(eventTable[k,"BLOCK"], eventTable[k,"id"])
    }
    
    # Compute epidemic variables x
    if ((nf <- length(f)) > 0L)
    {
        # check names(f) and initialize x-columns
        if (any(names(f) %in% names(data))) {
            warning("some 'names(f)' already existed in 'names(data)' ",
                    "and have been replaced")
        }
        data[names(f)] <- 0
        
        # Compute distance matrix
        firstDataBlock <- data[beginBlock[1L]:endBlock[1L],]
        coords <- firstDataBlock[coords.cols]
        rownames(coords) <- firstDataBlock[["id"]]
        distmat <- as.matrix(dist(coords, method = "euclidean"))
        diag(distmat) <- Inf   # no influence on yourself
        
        # Compute sum of distances over infectious individuals
        infectiousIDs <- firstDataBlock[firstDataBlock[["atRiskY"]] == 0, "id"]
        for(i in seq_along(beginBlock)) {
            blockidx <- beginBlock[i]:endBlock[i]
            blockdata <- data[blockidx,]
            blockIDs <- blockdata[["id"]]
            if (length(infectiousIDs) > 0L) {
                u <- distmat[as.character(blockIDs),as.character(infectiousIDs),
                             drop = FALSE]
                # numerical indices would be better, but be careful with factors
                data[blockidx,names(f)] <- sapply(f, function(B) rowSums(B(u)))
            }
            recoveredID <- blockIDs[blockdata[["Revent"]] == 1]
            infectedID <- blockIDs[blockdata[["event"]] == 1]
            if (length(recoveredID) > 0L) {
                infectiousIDs <- infectiousIDs[infectiousIDs != recoveredID]
            } else if (length(infectedID) > 0L) {
                infectiousIDs[length(infectiousIDs)+1L] <- infectedID
            }
        }
    }

    attr(data, "eventTimes") <- sort(eventTimes)
    attr(data, "timeRange") <- c(histIntervals[1L,1L],histIntervals[nBlocks,2L])
    attr(data, "coords.cols") <- coords.cols
    # must include this info because externally of this function
    # we don't know how many coords.cols (dimensions) we have
    attr(data, "f") <- f
    class(data) <- c("epidata", "data.frame")
    return(data)
}


################################################################################
# EXTRACTION OPERATOR FOR 'EPIDATA' OBJECTS
# Indexing with "[" would be possible (inheriting from data.frame).
# But using any column index would remove attributes (row indexes would not).
# Thus, we define an own method to retain and adjust the attributes when 
# selecting a subset of blocks of the 'epidata'.
# Selecting a subset of columns will remove class "epidata" (resulting in a
# simple data.frame)
################################################################################

"[.epidata" <- function(x, i, j, drop)
{
    # use data.frame method first
    xx <- NextMethod("[")
    # then return its result as pure data.frame or assure valid 'epidata'
    
    # if a subset of columns has been selected and attributes have been removed
    if (NCOL(xx) != ncol(x) || any(names(xx) != names(x))) {
        if (inherits(xx, "data.frame")) { # xx could be a vector
            class(xx) <- "data.frame"  # remove class 'epidata'
        }
        message("Note: converted class \"epidata\" to simple \"", class(xx),
                "\"")
        return(xx)
    }
    # else there is no effective column selection (e.g. j=TRUE)
    
    if (nrow(xx) == 0) {
        message("Note: no rows selected, dropped class \"epidata\"")
        class(xx) <- "data.frame"
        return(xx[TRUE])   # removes attributes
    }
    
    invalidEpidata <- FALSE
    blocksizesx <- table(x[["BLOCK"]])
    blocksizesxx <- table(xx[["BLOCK"]])
    blocksOK <- identical(c(blocksizesxx), c(blocksizesx[names(blocksizesxx)]))
    if (is.numeric(i) && any(diff(na.omit(i)) < 0)) {
        # epidata should remain ordered by time
        warning("dropped class \"epidata\": reordering rows is not permitted")
        invalidEpidata <- TRUE
    } else if (!blocksOK) {
        # blocks should not be cut, epidemic covariates might become invalid
        warning("dropped class \"epidata\": subsetting blocks not allowed")
        invalidEpidata <- TRUE
    } else if (any(diff(as.numeric(names(blocksizesxx))) != 1)) {
        # blocks can only be selected consecutively
        warning("dropped class \"epidata\": ",
                "only consecutive blocks may be selected")
        invalidEpidata <- TRUE
    }
    
    if (invalidEpidata) {
        class(xx) <- "data.frame"
        xx[TRUE] # removes attributes
    } else {
#         # adjust block index so that it starts at 1
#         firstBlockNumber <- as.numeric(names(blocksizesxx)[1])
#         if (firstBlockNumber > 1) {
#             xx[["BLOCK"]] <- xx[["BLOCK"]] - (firstBlockNumber-1)
#         }
        # Restore or adjust attributes
        tmin <- xx[["start"]][1]
        tmax <- xx[["stop"]][nrow(xx)]
        oldEventTimes <- attr(x, "eventTimes")
        attr(xx, "eventTimes") <-
            if (blocksOK) {
                oldEventTimes[oldEventTimes > tmin & oldEventTimes <= tmax]
            } else {
                xx[["stop"]][xx[["event"]] == 1]
            }
        attr(xx, "timeRange") <- c(tmin, tmax)
        attr(xx, "coords.cols") <- attr(x, "coords.cols")
        attr(xx, "f") <- attr(x, "f")
        xx
    }
}


################################################################################
# INSERT BLOCKS FOR EXTRA STOP TIMES IN 'EPIDATA' OBJECTS
################################################################################

intersperse <- function (epidata, stoptimes, verbose = FALSE)
{
    # Check arguments
    if (!inherits(epidata, "epidata")) {
        stop("'epidata' must inherit from class \"epidata\"")
    }
    if (!is.vector(stoptimes, mode = "numeric")) {
        stop("'stoptimes' must be a numeric vector")
    }
    
    # Identify new 'stoptimes'
    sortedEpiStop <- sort(unique(epidata$stop))
    extraStoptimes <- stoptimes[! stoptimes %in% sortedEpiStop]
    
    # Return original 'epidata' if nothing to do
    if (length(extraStoptimes) == 0) {
#         message("nothing done: no new stop times")
        return(epidata)
    }
    
#    # Retain attributes of 'epidata'
#    .attributes <- attributes(epidata)
#    .attributes <- .attributes[match(c("eventTimes", "timeRange",
#        "coords.cols", "f", "config", "call", "terms"), names(.attributes),
#        nomatch = 0)]

    # Check new 'stoptimes'
    timeRange <- attr(epidata, "timeRange")
    inside <- extraStoptimes > timeRange[1] & extraStoptimes < timeRange[2]
    if (any(!inside)) {
        extraStoptimes <- extraStoptimes[inside]
        warning("ignored extra 'stoptimes' outside the observation period")
    }
    
    # Impute blocks for extraStoptimes
    oldclass <- class(epidata)
    class(epidata) <- "data.frame" # Avoid use of [.epidata (not necessary here)
    blocksize <- sum(epidata$BLOCK == 1)
    nInsert <- length(extraStoptimes)
    lastRow <- nrow(epidata)
    epidata <- rbind(epidata,
                     epidata[rep.int(NA_integer_, nInsert * blocksize),],
                     deparse.level = 0) # add NA rows, to be replaced below
    if (verbose) pb <- txtProgressBar(min=0, max=nInsert, initial=0, style=3)
    for(i in seq_len(nInsert)) {
      extraStop <- extraStoptimes[i]
      nextStoptime <- sortedEpiStop[match(TRUE, sortedEpiStop > extraStop)]
      # Find the block (row indexes) into which the extraStop falls
      rowsMatchedBlock <- which(epidata$stop == nextStoptime)
      # Split this block up into 2 parts
      # later part equals original block with start time = extraStop
      newBlock <- epidata[rowsMatchedBlock,]
      newBlock$start <- extraStop
      # earlier part has stop time = extraStop and no events at this time point
      epidata[rowsMatchedBlock, "stop"] <- extraStop
      epidata[rowsMatchedBlock, "event"] <- 0
      epidata[rowsMatchedBlock, "Revent"] <- 0
      # write the new block to epidata (reorder rows later)
      epidata[lastRow + seq_along(rowsMatchedBlock),] <- newBlock
      lastRow <- lastRow + length(rowsMatchedBlock)
      if (verbose) setTxtProgressBar(pb, i)
    }
    if (verbose) close(pb)
    
    # Adjust BLOCK column
    sortedEpiStop <- sort(c(sortedEpiStop, extraStoptimes))
    epidata$BLOCK <- match(epidata$stop, sortedEpiStop)
    
    # Reorder rows by time and id
    epidata <- epidata[order(epidata$BLOCK, epidata$id), ]
    row.names(epidata) <- NULL
    class(epidata) <- oldclass
    
    return(epidata)
}


################################################################################
# SUMMARY FUNCTION FOR EPIDATA OBJECTS
# the epidemic is summarized by the following returned components:
# - type: one of "SIR", "SI", "SIRS", "SIS"
# - size: number of initially susceptible individuals, which became infected
# - initiallyInfected: vector (factor) of initially infected individuals
# - neverInfected: vector (factor) of never (during the observation period)
#                  infected individuals
# - coordinates: matrix with the coordinates of the individuals (rownames=id's)
# - byID: data.frame with time points of events by id (columns time.I, time.R
#         and optionally time.S)
# - counters: data.frame representing the evolution of the epidemic
################################################################################

summary.epidata <- function (object, ...)
{
    class(object) <- "data.frame"  # avoid use of [.epidata (not necessary here)
    
    # extract coordinates and initially infected individuals
    idlevels <- levels(object[["id"]])
    N <- length(idlevels)
    firstDataBlock <- object[object$BLOCK == min(object$BLOCK),]
    coordinates <- as.matrix(firstDataBlock[attr(object, "coords.cols")])
    rownames(coordinates) <- as.character(firstDataBlock[["id"]])
    initiallyInfected <- firstDataBlock$id[firstDataBlock$atRiskY == 0]
    m <- length(initiallyInfected)
    n <- N - m
    
    ### summary 1: event table with columns id, time and type (of event, S/I/R)
    # Extract time points of the S events for each id
    StimesID <- by(object[c("atRiskY", "stop")], object["id"],
                   function(x) {
                       SeventIdx <- which(diff(x[["atRiskY"]]) == 1)
                       x[["stop"]][SeventIdx]
                   }, simplify=TRUE)
    names(StimesID) <- paste0(names(StimesID), ":")
    StimesVec <- c(unlist(StimesID, use.names = TRUE)) # c() if by() returned an array
    .Sids <- sub("(.+):.*", "\\1", names(StimesVec))
    Stimes <- data.frame(id = factor(.Sids, levels = idlevels),
                         stop = StimesVec, type = rep("S", length(StimesVec)),
                         row.names = NULL, check.names = FALSE,
                         stringsAsFactors = FALSE)
    # Extract time points of the I and R events for each id
    Itimes <- object[object$event == 1, c("id", "stop")]
    Itimes[["type"]] <- rep("I", nrow(Itimes))
    Rtimes <- object[object$Revent == 1, c("id", "stop")]
    Rtimes[["type"]] <- rep("R", nrow(Rtimes))
    
    # Combine the three event types into one data.frame
    eventTable <- rbind(Rtimes, Stimes, Itimes)
      # need this order for the counters below in the case of SIS:
      # pseudo-R-event occures infinitesimally before S
    names(eventTable)[2L] <- "time"
    eventTable <- eventTable[order(eventTable[["id"]], eventTable[["time"]]), ]
    eventTable[["type"]] <- factor(eventTable[["type"]], levels=c("S","I","R"))
    rownames(eventTable) <- NULL
    
    ### summary 2: type and size of the epidemic
    resusceptibility <- length(StimesVec) > 0
    epitype <-
        if (resusceptibility) {
            Rtimes_notLast <- Rtimes[-which.max(Rtimes[,2]),]
            onlyPseudoR <- length(setdiff(Rtimes_notLast[,2], Stimes[,2])) == 0
            if (onlyPseudoR) "SIS" else "SIRS"
        } else {
            if (nrow(Rtimes) > 0) "SIR" else "SI"
        }
    isEverInfected <- idlevels %in% initiallyInfected |
        idlevels %in% unique(eventTable$id[eventTable$type == "I"])
    isNeverInfected <- !isEverInfected
    size <- n - sum(isNeverInfected)
#     everInfected <- factor(idlevels[isEverInfected], levels = idlevels)
    neverInfected <- factor(idlevels[isNeverInfected], levels = idlevels)
    
    ### summary 3: eventTable by id in wide form
    byID_everInfected <-
        if (nrow(eventTable) == 0) {
            data.frame(id = factor(character(0), levels = idlevels),
                       time.I = numeric(0), row.names = NULL,
                       check.names = FALSE, stringsAsFactors = FALSE)
        } else if (!resusceptibility) {
            .res <- reshape(eventTable, direction = "wide", timevar = "type",
                           idvar = "id")
            attr(.res, "reshapeWide") <- NULL
            .res
        } else {
            rowsPerId <- table(eventTable[["id"]])
            modulo3 <- rowsPerId %% 3
            rest1 <- modulo3 == 1
            rest12 <- modulo3 >= 1
            missingR <-
                data.frame(id = names(rowsPerId)[rest1],
                           time = rep(NA_real_, sum(rest1)),
                           type = rep("R", sum(rest1)), row.names = NULL,
                           check.names = FALSE, stringsAsFactors = FALSE)
            missingS <- 
                data.frame(id = names(rowsPerId)[rest12],
                           time = rep(NA_real_, sum(rest12)),
                           type = rep("S", sum(rest12)), row.names = NULL,
                           check.names = FALSE, stringsAsFactors = FALSE)
            eventTable3 <- rbind(eventTable, missingR, missingS)
            eventTable3 <- eventTable3[order(eventTable3[["id"]]),]
            .res <- data.frame(
                eventTable3[eventTable3$type == "I", c("id", "time")],
                eventTable3[eventTable3$type == "R", "time", drop = FALSE],
                eventTable3[eventTable3$type == "S", "time", drop = FALSE],
                row.names = NULL, check.names = FALSE, stringsAsFactors = FALSE
            )
            names(.res) <- c("id", paste("time", c("I", "R", "S"), sep="."))
            .res
        }
    byID_neverInfected <- data.frame(id = neverInfected,
        time.I = rep(NA_real_, n-size), time.R = rep(NA_real_, n-size),
        time.S = rep(NA_real_, n-size), row.names = NULL, check.names = FALSE)
    byID_all <- rbind(byID_everInfected,
                      byID_neverInfected[seq_along(byID_everInfected)])
    byID <- byID_all[order(byID_all[["id"]]),]
    rownames(byID) <- NULL
    
    ### summary 4: upgrade eventTable with
    ###            evolution of numbers of susceptibles, infectious and removed
    counters <- eventTable[order(eventTable[["time"]]),c("time", "type", "id")]
    init <- data.frame(time = attr(object, "timeRange")[1L],
                       type = NA_character_, id = NA_character_,
                       nSusceptible = n, nInfectious = m, nRemoved = 0L)
    cumulatedReSusceptibility <- cumsum(counters[["type"]] == "S")
    cumulatedInfections <- cumsum(counters[["type"]] == "I")
    cumulatedRemovals <- cumsum(counters[["type"]] == "R")
    counters[["nSusceptible"]] <-
        init[["nSusceptible"]] - cumulatedInfections + cumulatedReSusceptibility
    counters[["nInfectious"]] <-
        init[["nInfectious"]]  + cumulatedInfections - cumulatedRemovals
    counters[["nRemoved"]] <-
        init[["nRemoved"]]     + cumulatedRemovals   - cumulatedReSusceptibility
    counters <- rbind(init, counters)
    rownames(counters) <- NULL

    ### return the components in a list
    res <- list(type = epitype, size = n - sum(isNeverInfected),
        initiallyInfected = initiallyInfected, neverInfected = neverInfected,
        coordinates = coordinates, byID = byID, counters = counters)
    class(res) <- "summary.epidata"
    attr(res, "eventTimes") <- attr(object, "eventTimes")
    attr(res, "timeRange") <- attr(object, "timeRange")
    res
}



################################################################################
# PRINT METHOD FOR 'EPIDATA' OBJECTS
################################################################################

print.epidata <- function (x, ...)
{
    cat("\nHistory of an epidemic\n")
    cat("Number of individuals:", nlevels(x[["id"]]), "\n")
    cat("Time range:", paste(attr(x, "timeRange"), collapse = " -- "), "\n")
    cat("Number of infections:", length(attr(x, "eventTimes")), "\n\n")
    print.data.frame(x, ...)
    cat("\n")
    invisible(x)
}


################################################################################
# PRINT METHOD FOR THE SUMMARY OF 'EPIDATA' OBJECTS
################################################################################

print.summary.epidata <- function(x, ...)
{
    cat("\nAN", x$type, "EPIDEMIC\n")
    cat("  Time range:", paste(attr(x, "timeRange"), collapse = " -- "), "\n")
    cat("  Number of individuals:", nlevels(x$initiallyInfected), "\n")
    cat(" ", length(x$initiallyInfected), "initially infected individuals")
    if (length(x$initiallyInfected) > 0) {
        cat(":\n    ")
        str(as.character(x$initiallyInfected), give.head = FALSE, vec.len = 100,
            strict.width = "wrap", indent.str = "  ")
    } else cat("\n")
    cat(" ", length(x$neverInfected), "never infected individuals")
    if (length(x$neverInfected) > 0) {
        cat(":\n    ")
        str(as.character(x$neverInfected), give.head = FALSE, vec.len = 100,
            strict.width = "wrap", indent.str = "  ")
    } else cat("\n")
    cat("  Size of the epidemic:", x$size, "\n")
    if (x$type %in% c("SIRS", "SIS")) {
        cat("  Number of infections:", length(attr(x, "eventTimes")), "\n")
    }
    dimc <- dim(x$counters)
    cat("\n$ counters ('data.frame',", dimc[1L], "x", dimc[2L], "):",
        "evolution of the epidemic:\n")
    counters2print <- if (dimc[1] > 6L) {
            tmp <- format.data.frame(x$counters[c(1:3,1,dimc[1]-(1:0)),],
                                     na.encode = FALSE)
            tmp[4,] <- c("[....]", "", "", "", "", "")
            rownames(tmp)[4] <- ""
            as.matrix(tmp)
        } else { x$counters }
    print(counters2print, quote = FALSE, right = TRUE, na.print = "")
    cat("\n")
    invisible(x)
}


################################################################################
# PLOT METHOD FOR "summary.epidata" OR "epidata" OBJECTS
# shows the evolution of the numbers of susceptible, infectious and recovered
# individuals.
################################################################################

plot.summary.epidata <- function (x,
    lty = c(2,1,3), lwd = 2, col = 1, col.hor = col, col.vert = col,
    xlab = "Time", ylab = "Number of individuals", xlim = NULL, ylim = NULL,
    legend.opts = list(), do.axis4 = NULL, panel.first = grid(),
    rug.opts = list(), which.rug = c("infections", "removals",
    "susceptibility", "all"), ...)
{
    counters <- x[["counters"]]
    type <- x[["type"]]

    n <- counters[1L,"nSusceptible"]
    m <- counters[1L,"nInfectious"]
    N <- n + m
    times <- counters[-1L,"time"]
    if (missing(lty)) {
        lty <- c(2, 1, 3 * (type %in% c("SIR","SIRS")))
    }
    recycle3 <- function (xnam)
        assign(xnam, rep(get(xnam), length.out = 3), inherits = TRUE)
    for(varname in c("lty", "lwd", "col", "col.hor", "col.vert"))
        recycle3(varname)
    if (is.null(xlim)) {
        xlim <- attr(x, "timeRange")
        if (xlim[2] == Inf) xlim[2] <- times[length(times)]
    }
    if (is.null(ylim))
        ylim <- c(0, max(
            (lty[1] > 0) * {if (type %in% c("SIRS", "SIS")) N else n},
            (lty[2] > 0) * max(counters$nInfectious),
            (lty[3] > 0) * max(counters$nRemoved)
            ))

    # basic plotting frame
    plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab,
         panel.first = panel.first, ...)
    abline(h = c(0, N), col = "grey")

    # for real xlim in lines.stepfun (see 'dr' adjustment in plot.stepfun code)
    fakexlim <- c(1,2) * (xlim[2] + 2*xlim[1])/3 - c(0,xlim[1])
    # this isn't nice, a user argument 'dr' in plot.stepfun would be appreciated
    
    # add #Susceptibles
    if (all(counters$nSusceptible == n)) {
        lines(x = xlim, y = c(n,n),
              lty = lty[1], lwd = lwd[1], col = col.hor[1], ...)
    } else {
        lines(stepfun(times, counters$nSusceptible), xlim = fakexlim,
              lty = lty[1], lwd = lwd[1], col.hor = col.hor[1],
              col.vert = col.vert[1], do.points = FALSE, ...)
    }

    # add #Infected
    if (all(counters$nInfectious == m)) {
        lines(x = xlim, y = c(m,m),
              lty = lty[2], lwd = lwd[2], col = col.hor[2], ...)
    } else {
        lines(stepfun(times, counters$nInfectious),
              xlim = fakexlim, lty = lty[2], lwd = lwd[2], col.hor = col.hor[2],
              col.vert = col.vert[2], do.points = FALSE, ...)
    }

    # add #Removed
    if (all(counters$nRemoved == 0)) {
        lines(x = xlim, y = c(0,0),
              lty = lty[3], lwd = lwd[3], col = col.hor[3], ...)
    } else {
        lines(stepfun(times, counters$nRemoved),
              xlim = fakexlim, lty = lty[3], lwd = lwd[3], col.hor = col.hor[3],
              col.vert = col.vert[3], do.points = FALSE, ...)
    }

    # add special annotations
    if (is.null(do.axis4)) do.axis4 <- type == "SIR"
    if (do.axis4) {
        finalvalues <- counters[nrow(counters), c("nSusceptible", "nRemoved")]
        axis(4, at = finalvalues[lty[c(1,3)] > 0], font = 2, ...)
    }
    if (is.list(rug.opts)) {
        if (is.null(rug.opts$ticksize)) rug.opts$ticksize <- 0.02
        if (is.null(rug.opts$quiet)) rug.opts$quiet <- TRUE
        which.rug <- match.arg(which.rug)
        if (is.null(rug.opts$col)) rug.opts$col <-
            switch(which.rug, all = 1, infections = col.hor[2],
                   removals = col.hor[3], susceptibility = col.hor[1])
        rugLocations <- switch(which.rug,
            all = times, infections = attr(x, "eventTimes"),
            removals = counters$time[counters$type == "R"],
            susceptibility = counters$time[counters$type == "S"]
        )
        if (length(rugLocations) > 0) {
            do.call(rug, c(list(x = rugLocations), rug.opts))
        }
    }
    if (is.list(legend.opts)) {
        if (is.null(legend.opts[["x",exact=TRUE]]))
            legend.opts$x <- "topright"
        if (is.null(legend.opts$legend))
            legend.opts$legend <- c("susceptible", "infectious", "removed")
        if (any(lty == 0) && length(legend.opts$legend) == 3)
            legend.opts$legend <- legend.opts$legend[lty > 0]
        if (is.null(legend.opts$lty)) legend.opts$lty <- lty[lty > 0]
        if (is.null(legend.opts$lwd)) legend.opts$lwd <- lwd[lty > 0]
        if (is.null(legend.opts$col)) legend.opts$col <- col.hor[lty > 0]
        if (is.null(legend.opts$bty)) legend.opts$bty <- "n"
        do.call(legend, legend.opts)
    }
    invisible(as.matrix(
        counters[c("time", "nSusceptible", "nInfectious", "nRemoved")]
    ))
}

plot.epidata <- function(x, ...)
{
    sx <- summary(x)
    plot(sx, ...)
}


################################################################################
# ANIMATION OF EPIDEMICS
# spatio-temporal animation (for 1/2d-coordinates only)
# two types:
#   - sequential plots regardless of time between events (i.e. only ordering)
#   - chronological animation with timer
################################################################################

animate.summary.epidata <- function (object,
    main = "An animation of the epidemic",
    pch = 19, col = c(3, 2, gray(0.6)), time.spacing = NULL,
    sleep = quote(5/.nTimes), legend.opts = list(), timer.opts = list(),
    end = NULL, generate.snapshots = NULL, ...)
{
    counters <- object[["counters"]]
    # remove pseudo-R-events, which come before S-event
    directSevents <- which(duplicated(counters[["time"]]))
    counters_noPseudoR <- if (length(directSevents)) {
            counters[-(directSevents-1), ]
        } else {
            counters
        }
    # remove initial row and keep essential columns
    eventTable <- counters_noPseudoR[-1, c("time", "type", "id")]
    eventTable[["type"]] <- unclass(eventTable[["type"]])  # get integer codes
    .nTimes <- nrow(eventTable)
    
    # extract initial individual information (id, at-risk, coordinates)
    coords <- object[["coordinates"]]
    d <- ncol(coords)
    if (d > 2L) {
        stop("spatial plotting in more than two dimensions is not implemented")
    } else if (d == 1L) {
        coords <- cbind(coords, 0)
    } else if (d == 0L) {
        stop ("'object' does not contain any defined coordinates")
    }
    
    # plot the initial state
    pch <- rep(pch, length.out = 3)
    col <- rep(col, length.out = 3)
    isInitiallyInfected <- rownames(coords) %in% object[["initiallyInfected"]]
    plot(coords, pch = ifelse(isInitiallyInfected, pch[2L], pch[1L]), 
                 col = ifelse(isInitiallyInfected, col[2L], col[1L]),
                 main = main, ...)
    if (is.list(legend.opts)) {
        if (is.null(legend.opts[["x",exact=TRUE]]))
            legend.opts$x <- "topright"
        if (is.null(legend.opts$legend))
            legend.opts$legend <- c("susceptible", "infectious", "removed")
        if (is.null(legend.opts$col)) legend.opts$col <- col
        if (is.null(legend.opts$pch)) legend.opts$pch <- pch
        do.call(legend, legend.opts)
    }
    
    # animate the epidemic by iteratively re-drawing points at the coordinates
    sleep <- eval(sleep)
    if (is.null(time.spacing)) { # plot events sequentially
        for(i in seq_len(.nTimes)) {
            Sys.sleep(sleep)
            tmp <- eventTable[i,]  # c(time, type, id)
            points(coords[as.character(tmp[["id"]]),,drop=FALSE],
                   pch = pch[tmp[["type"]]], col = col[tmp[["type"]]])
        }
    } else { # plot events chronologically
        if (is.null(end))
            end <- eventTable[.nTimes, "time"] + time.spacing
        timeGrid <- seq(from = time.spacing, to = end, by = time.spacing)
        timeWidth <- nchar(timeGrid[length(timeGrid)])
        timeDigits <- nchar(strsplit(as.character(time.spacing), ".",
            fixed = TRUE)[[1L]][2L])
        form <- paste("%", timeWidth, ".", timeDigits, "f", sep = "")
        if (is.list(timer.opts)) {
            if (is.null(timer.opts[["x",exact=TRUE]]))
                timer.opts$x <- "bottomright"
            if (is.null(timer.opts$title))   timer.opts$title <- "time"
            if (is.null(timer.opts$box.lty)) timer.opts$box.lty <- 0
            if (is.null(timer.opts$adj))     timer.opts$adj <- c(0.5,0.5)
            if (is.null(timer.opts$inset))   timer.opts$inset <- 0.01
            if (is.null(timer.opts$bg))      timer.opts$bg <- "white"
            do.call(legend, c(list(legend = sprintf(form, 0)), timer.opts))
        }
        oldtp <- tp <- attr(object, "timeRange")[1L]
        i <- 1L                   # to be used in the file argument in dev.print
        if (is.vector(generate.snapshots, mode="character") &&
            length(generate.snapshots) == 1L && require("animation")) {
            img.name <- generate.snapshots
            ani.dev <- animation::ani.options("ani.dev")
            if (is.character(ani.dev)) ani.dev <- get(ani.dev)
            imgdir <- file.path(animation::ani.options("outdir"),
                                animation::ani.options("imgdir"))
            imgtype <- animation::ani.options("ani.type")
            generate.snapshots <- list(
                device = ani.dev,
                file = quote(file.path(imgdir, paste0(img.name,i,".",imgtype))),
                width = animation::ani.options("ani.width"),
                height = animation::ani.options("ani.height")
            )
        }
        if (is.list(generate.snapshots)) {
            do.call(dev.print, generate.snapshots)
        }
        for(i in 1L+seq_along(timeGrid)) {
            tp <- timeGrid[i-1L]
            Sys.sleep(sleep)
            timeIndex <- which(eventTable[["time"]] > oldtp & eventTable[["time"]] <= tp)
            if (length(timeIndex) > 0L) {
                tmp <- eventTable[timeIndex,]  # c(time, type, id)
                points(coords[as.character(tmp[["id"]]),,drop=FALSE],
                       pch = pch[tmp[["type"]]], col = col[tmp[["type"]]])
            }
            if (is.list(timer.opts)) {
                do.call(legend, c(list(legend = sprintf(form,tp)), timer.opts))
            }
            oldtp <- tp
            if (is.list(generate.snapshots)) {
                do.call(dev.print, generate.snapshots)
            }
        }
    }
    invisible(NULL)
}

animate.epidata <- function (object, ...)
{
    s <- summary(object)
    animate(s, ...)
}



################################################################################
# PLOT THE STATE CHANGES OF INDIVIDUALS IN AN EPIDEMIC ('EPIDATA' OBJECT)
# ... will be passed to the plot function (stepfun or curve),
# e.g. add, xlim, ylim, main, xlab, ylab, ...
################################################################################

stateplot <- function(x, id, ...)
{
    sx <- getSummary(x, class = "epidata")

    .id <- as.character(id)
    if (length(.id) != 1) {
        stop ("'id' must have length 1")
    }
    initiallyInfected <- sx[["initiallyInfected"]]
    if (! .id %in% levels(initiallyInfected)) {
        stop ("invalid 'id', does not exist in 'x'")
    }
    isInitiallyInfected <- .id %in% initiallyInfected
    
    counters <- sx[["counters"]]
    states <- levels(counters[["type"]])
    
    path <- counters[which(counters$id == .id), c("time", "type")]
    # remove pseudo-R-events, which come before S-event
    directSevents <- which(duplicated(path[["time"]]))
    path_noPseudoR <- if (length(directSevents)) {
            path[-(directSevents-1), ]
        } else {
            path
        }
    
    pathfunc <-
        if (nrow(path_noPseudoR) > 0) {
            stepfun(
                x = path_noPseudoR[["time"]],
                y = c(1+isInitiallyInfected, unclass(path_noPseudoR[["type"]])),
                right = FALSE
            )
        } else {
            function(t) rep(1+isInitiallyInfected, length(t))
        }
    
    # plot it
    dotargs <- list(...)
    nms <- names(dotargs)
    if(! "xlab" %in% nms) dotargs$xlab <- "time"
    if(! "ylab" %in% nms) dotargs$ylab <- "state"
    if(! "main" %in% nms) dotargs$main <- ""
    if(! "xlim" %in% nms) dotargs$xlim <- attr(sx, "timeRange")
    if(! "xaxs" %in% nms) dotargs$xaxs <- "i"
    if(! "do.points" %in% nms && inherits(pathfunc, "stepfun")) {
        dotargs$do.points <- FALSE
    }
    do.call("plot", args = c(list(x = pathfunc, yaxt = "n"), dotargs))
    axis(2, at = seq_along(states), labels = states)

    invisible(pathfunc)
}
