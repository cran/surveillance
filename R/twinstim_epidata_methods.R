################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### S3-methods for "epidataCS" data objects, which represent
### CONTINUOUS SPATIO-temporal infectious disease case data
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 666 $
### $Date: 2013-11-08 15:45:36 +0100 (Fre, 08 Nov 2013) $
################################################################################


### Number of events

nobs.epidataCS <- function (object, ...) length(object$events)


### UPDATE eps.s, eps.t, qmatrix OR nCircle2Poly IN AN EXISTING epidataCS OBJECT

# all arguments but 'object' are optional, the ... argument is unused
update.epidataCS <- function (object, eps.t, eps.s, qmatrix, nCircle2Poly, ...)
{
    nEvents <- nobs(object)

    # Check and update eps.t
    if (!missing(eps.t)) {
        stopifnot(is.numeric(eps.t), eps.t > 0)
        object$events$eps.t <- eps.t
    }

    # Initialise indicator of which influenceRegions to update
    ir2update <- logical(nEvents)   # all FALSE

    # Check and update eps.s
    if (!missing(eps.s)) {
        stopifnot(is.numeric(eps.s), eps.s > 0)
        oldeps.s <- object$events$eps.s
        object$events$eps.s <- eps.s
        ir2update <- oldeps.s != object$events$eps.s
    }

    # Check nCircle2Poly
    nCircle2Poly <- if (missing(nCircle2Poly)) {
        attr(object$events$.influenceRegion, "nCircle2Poly")
    } else {
        stopifnot(isScalar(nCircle2Poly))
        ir2update <- rep.int(TRUE, nEvents)
        as.integer(nCircle2Poly)
    }

    # Update influenceRegions of events
    if (any(ir2update)) {
        object$events$.influenceRegion[ir2update] <-
            .influenceRegions(object$events[ir2update,], object$W, nCircle2Poly)
        attr(object$events$.influenceRegion, "nCircle2Poly") <- nCircle2Poly
    }

    # Check qmatrix
    if (!missing(qmatrix)) object$qmatrix <- checkQ(qmatrix, levels(object$events$type))

    #hoehle @ 16 Apr 2011 - bug fix. .obsInfLength was not handled
    # Update length of infection time, i.e. length = min(T-time, eps.t)
    if (!missing(eps.t)) {
      timeRange <- with(object$stgrid, c(start[1], stop[length(stop)]))
      object$events$.obsInfLength <- with(object$events@data, pmin(timeRange[2]-time, eps.t))
    }

    # Update .sources
    if (!missing(eps.t) || !missing(eps.s) || !missing(qmatrix)) {
        object$events$.sources <- determineSources.epidataCS(object)
    }

    # Done update.
    return(object)
}



### subsetting epidataCS, i.e. select only part of the events,
### but retain stgrid and W. If any event types disappear due to subsetting,
### these types will be dropped from the factor levels and from qmatrix

"[.epidataCS" <- function (x, i, j, drop = FALSE)
{
    ## Store nCircle2Poly attribute of x$events$.influenceRegion since this will
    ## be dropped when subsetting
    nCircle2Poly <- attr(x$events$.influenceRegion, "nCircle2Poly")

    ## apply [,SpatialPointsDataFrame-method
    cl <- sys.call()
    cl[[1]] <- as.name("[")
    cl[[2]] <- substitute(x$events)
    x$events <- eval(cl, envir=parent.frame())

    ## restore nCircle2Poly attribute
    attr(x$events$.influenceRegion, "nCircle2Poly") <- nCircle2Poly
    
    ## assure valid epidataCS after subsetting
    if (!missing(j)) {                # only epidemic covariates may be selected
        BLOCKstartEndemicVars <- setdiff(names(x$stgrid),
            setdiff(obligColsNames_stgrid,"start"))
        if (!all(obligColsNames_events %in% names(x$events)) ||
            !all(BLOCKstartEndemicVars %in% names(x$events))) {
            stop("only epidemic covariates may be removed from 'events'")
        }
    }
    if (!missing(i)) {
        ## update .sources
        x$events$.sources <- determineSources.epidataCS(x)
        ## update types and qmatrix (a type could have disappeared)
        x$events$type <- x$events$type[drop=TRUE]
        typeNames <- levels(x$events$type)
        if (!identical(rownames(x$qmatrix), typeNames)) {
            message("Note: dropped type(s) ",
                    paste0("\"", setdiff(rownames(x$qmatrix), typeNames), "\"",
                           collapse = ", "))
            typesIdx <- match(typeNames, rownames(x$qmatrix))
            x$qmatrix <- x$qmatrix[typesIdx, typesIdx, drop = FALSE]
        }
    }

    ## Done
    return(x)
}


## The subset method for epidataCS-objects is adapted from
## base::subset.data.frame (authored by Peter 
## Dalgaard and Brian Ripley, Copyright (C) 1995-2012
## The R Core Team) with slight modifications only
## (we just replace 'x' by 'x$events@data' for evaluation of subset and select)

subset.epidataCS <- function (x, subset, select, drop = FALSE, ...)
{
    if (missing(subset)) 
        r <- TRUE
    else {
        e <- substitute(subset)
        r <- eval(e, x$events@data, parent.frame()) # HERE IS A MOD
        if (!is.logical(r)) stop("'subset' must evaluate to logical")
        r <- r & !is.na(r)
    }
    if (missing(select)) 
        vars <- TRUE
    else {
        nl <- as.list(seq_along(x$events@data)) # HERE IS A MOD
        names(nl) <- names(x$events@data)       # HERE IS A MOD
        vars <- eval(substitute(select), nl, parent.frame())
    }
    x[r, vars, drop = drop]       # this calls the [.epidataCS-method from above
}


## Subset epidataCS object using head and tail methods (which use [.epidataCS)

head.epidataCS <- function (x, n = 6L, ...)
    head.matrix(x, n = n, ...)

tail.epidataCS <- function (x, n = 6L, ...)
{
    # ugly hack for tail.matrix because I don't want to register a
    # dim-method for class "epidataCS"
    nrow <- function (x) base::nrow(x$events)
    my.tail.matrix <- tail.matrix
    environment(my.tail.matrix) <- environment()
    ##<- such that the function uses my local nrow definition
    my.tail.matrix(x, n = n, addrownums=FALSE, ...)
}



### extract marks of the events (actually also including time and tile)

marks.epidataCS <- function (x, coords = TRUE, ...) {
    noEventMarks <- setdiff(reservedColsNames_events, "ID")
    endemicCovars <- setdiff(names(x$stgrid),
        c(reservedColsNames_stgrid, obligColsNames_stgrid))
    idxnonmarks <- match(c(noEventMarks, endemicCovars), names(x$events))
    if (coords) {          # use as.data.frame method for SpatialPointsDataFrame
        as.data.frame(x$events[-idxnonmarks])
    } else {                            # return marks without coordinates
        x$events@data[-idxnonmarks]
    }
}



### printing methods

print.epidataCS <- function (x, n = 6L, digits = getOption("digits"), ...)
{
    timeRange <- c(x$stgrid$start[1], x$stgrid$stop[nrow(x$stgrid)])
    bboxtxt <- paste(apply(bbox(x$W), 1,
        function (int) paste0("[", paste(format(int, trim=TRUE, digits=digits), collapse=", "), "]")
        ), collapse = " x ")
    nBlocks <- length(unique(x$stgrid$BLOCK))
    nTiles <- nlevels(x$stgrid$tile)
    typeNames <- levels(x$events$type)
    nEvents <- nobs(x)
    cat("\nHistory of an epidemic\n")
    cat("Observation period:", paste(format(timeRange, trim=TRUE, digits=digits), collapse = " -- "), "\n")
    cat("Observation window (bounding box):", bboxtxt, "\n")
    cat("Spatio-temporal grid (not shown):", nBlocks,
        ngettext(nBlocks, "time block,", "time blocks,"),
        nTiles, ngettext(nTiles, "tile", "tiles"), "\n")
    cat("Types of events:", paste0("'",typeNames,"'"), "\n")
    cat("Overall number of events:", nEvents, "\n\n")
    # 'print.SpatialPointsDataFrame' does not pass its "digits" argument on to
    # 'print.data.frame', hence the use of options()
    odigits <- options(digits=digits); on.exit(options(odigits))
    visibleCols <- grep("^\\..+", names(x$events@data), invert = TRUE)
    print(head.matrix(x$events[visibleCols], n = n), ...)
    if (n < nEvents) cat("[....]\n")
    cat("\n")
    invisible(x)
}



### SUMMARY
# the epidemic is summarized by the following returned components:
# timeRange, nEvents, eventTimes, eventCoords, nSources, as well as
# - tile/typetable: number of events per tile/type
# - counter: number of infective individuals as stepfun

summary.epidataCS <- function (object, ...)
{
    coords <- coordinates(object$events)
    times <- object$events$time
    nEvents <- length(times)
    timeRange <- with(object$stgrid, c(start[1], stop[length(stop)]))
    nBlocks <- length(unique(object$stgrid$BLOCK))
    tiles <- object$events$tile
    bbox <- bbox(object$W)
    tileTable <- c(table(tiles))
    types <- object$events$type
    nTypes <- nlevels(types)
    typeTable <- c(table(types))
    nSources <- sapply(object$events$.sources, length)
    eps <- object$events@data[c("eps.t", "eps.s")]
    eventMarks <- marks(object)

    removalTimes <- times + object$events$eps.t
    tps <- sort(unique(c(times, removalTimes[is.finite(removalTimes)])))
    nInfectious <- sapply(tps, function(t) sum(times <= t & removalTimes > t))
    counter <- stepfun(tps, c(0,nInfectious), right = TRUE)

    res <- list(timeRange = timeRange, bbox = bbox, nBlocks = nBlocks,
                nEvents = nEvents, nTypes = nTypes,
                eventTimes = times, eventCoords = coords, eventTypes = types,
                eventRanges = eps, eventMarks = eventMarks, 
                tileTable = tileTable, typeTable = typeTable,
                counter = counter, nSources = nSources)
    class(res) <- "summary.epidataCS"
    res
}

print.summary.epidataCS <- function (x, ...)
{
    bboxtxt <- paste(apply(x$bbox, 1,
        function (int) paste0("[", paste(format(int, trim=TRUE), collapse=", "), "]")
        ), collapse = " x ")
    cat("\n")
    cat("Observation period:", paste(format(x$timeRange, trim=TRUE), collapse = " -- "), "\n")
    cat("Observation window (bounding box):", bboxtxt, "\n")
    cat("Spatio-temporal grid (not shown):", x$nBlocks,
        ngettext(x$nBlocks, "time block,", "time blocks,"),
        length(x$tileTable), ngettext(length(x$tileTable), "tile", "tiles"), "\n")
    cat("Overall number of events:", x$nEvents,
        if (x$nTypes==1) "(single type)" else paste0("(",x$nTypes," types)"),
        "\n")
    
    ## if (x$nTypes > 1) {
    ##     cat(x$nTypes, "types of events:\n")
    ##     print(as.table(x$typeTable))
    ## }
    ## cat("\nTime points of events:\n")
    ## print(summary(x$eventTimes))
    ## cat("\nTiles of events:\n")
    ## print(as.table(x$tileTable))
    cat("\nSummary of the event marks:\n")
    print(summary(x$eventMarks))

    cat("Number of potential sources of transmission:\n")
    if (any(is.finite(unlist(x$eventRanges)))) {
        print(summary(x$nSources))
        cat("\nStep function of the number of infectives over time:\n")
        print(x["counter"])
    } else {
        cat("   monotonically increasing like the number of infectives\n",
            "   because of infinite ranges of interaction ('eps.t' and 'eps.s')\n", sep="")
    }
    cat("\n")
    
    invisible(x)
}



### animate
# spatio-temporal animation, two types:
# time.spacing=NULL: sequential plots regardless of time between events (i.e. only ordering)
# time.spacing=scalar: chronological animation with timer. if time.spacing = NA, then the time step is automatically determined such that ani.options("nmax") snapshots result.
# respects ani.options "interval" and "nmax"

animate.epidataCS <- function (object, interval = c(0,Inf), time.spacing = NULL,
    nmax = NULL, sleep = NULL, legend.opts = list(), timer.opts = list(),
    pch = 15:18, col.current = "red", col.I = "#C16E41", col.R = "#B3B3B3",
    col.influence = "#FEE0D2", main = NULL, verbose = TRUE, ...)
{
    stopifnot(is.numeric(interval), length(interval) == 2L)
    with.animation <- suppressWarnings(require("animation"))
    if (is.null(sleep)) {
        sleep <- if (with.animation) animation::ani.options("interval") else 0.1
        ## we cannot set this as default function argument, because we don't
        ## want to depend on package "animation" (surveillance only suggests it)
    }
    if (is.null(nmax)) {
        nmax <- if (with.animation) animation::ani.options("nmax") else Inf
    }
    s <- summary(object)
    removalTimes <- s$eventTimes + object$events$eps.t
    eventCoordsTypes <- cbind(s$eventCoords, type = s$eventTypes)
    pch <- rep_len(pch, s$nTypes)
    typeNames <- names(s$typeTable)
    multitype <- length(typeNames) > 1L

    # set default legend options
    doLegend <- if (is.list(legend.opts)) {
        if (is.null(legend.opts[["x"]])) legend.opts$x <- "topright"
        if (is.null(legend.opts$title))  legend.opts$title <-
            if (multitype) "type" else "state"
        if (is.null(legend.opts$legend)) { legend.opts$legend <-
            if (multitype) typeNames else c("infectious", if (!is.na(col.R)) "removed")
        }
        if (is.null(legend.opts$col)) { legend.opts$col <-
            if (multitype) col.current else c(col.I, if (!is.na(col.R)) col.R)
        }
        if (is.null(legend.opts$pch)) legend.opts$pch <- pch
        TRUE
    } else FALSE

    # set default timer options
    doTimer <- if (is.list(timer.opts)) {
        if (is.null(timer.opts[["x"]]))  timer.opts$x <- "bottomright"
        if (is.null(timer.opts$title))   timer.opts$title <- "time"
        if (is.null(timer.opts$box.lty)) timer.opts$box.lty <- 0
        if (is.null(timer.opts$adj))     timer.opts$adj <- c(0.5,0.5)
        if (is.null(timer.opts$inset))   timer.opts$inset <- 0.01
        if (is.null(timer.opts$bg))      timer.opts$bg <- "white"
        TRUE
    } else FALSE

    # determines multiplicity of rows of a numeric matrix
    # and returns unique rows with appended multiplicity column
    countunique <- function (mat) {
        count <- multiplicity(mat)
        unique(cbind(mat, count))
    }
    # wrapper for 'points' with specific 'cex' for multiplicity
    multpoints <- function (tableCoordsTypes, col) {
        tableMult <- countunique(tableCoordsTypes)
        points(tableMult[,1:2,drop=FALSE], pch = pch[tableMult[,"type"]],
               col = col, cex = sqrt(1.5*tableMult[,"count"]/pi) * par("cex"))
    }
    # functions returning if events are in status I or R at time t
    I <- function (t) s$eventTimes <= t & removalTimes >= t
    R <- function (t) removalTimes < t

    sequential <- is.null(time.spacing)  # plot observed infections sequentially
    if (!sequential) stopifnot(length(time.spacing) == 1L)
    timeGrid <- if (sequential) unique(s$eventTimes) else {
        start <- max(s$timeRange[1], interval[1])
        end <- min(interval[2], s$timeRange[2],
            max(removalTimes) + if (is.na(time.spacing)) 0 else time.spacing)
        if (is.na(time.spacing)) {
            if (!is.finite(nmax)) {
                stop("with 'time.spacing=NA', 'nmax' must be finite")
            }
            seq(from = start, to = end, length.out = nmax)
        } else {
            tps <- seq(from = start, to = end, by = time.spacing)
            if (length(tps) > nmax) {
                message("Generating only the first ",
                        sQuote(if (with.animation) "ani.options(\"nmax\")" else "nmax"),
                        " (=", nmax, ") snapshots")
                head(tps, nmax)
            } else tps
        }
    }
    .info <- format.info(timeGrid)
    timerformat <- paste0("%", .info[1], ".", .info[2], "f")

    # animate
    loopIndex <- if (!sequential) timeGrid else {
        idxs <- which(s$eventTimes >= interval[1] & s$eventTimes <= interval[2])
        if (length(idxs) > nmax) {
            message("Generating only the first ",
                    sQuote(if (with.animation) "ani.options(\"nmax\")" else "nmax"),
                    " (=", nmax, ") events")
            head(idxs, nmax)
        } else idxs
    }
    told <- -Inf
    if (verbose)
        pb <- txtProgressBar(min=0, max=max(loopIndex), initial=0, style=3)
    for(it in loopIndex) {
        t <- if (sequential) s$eventTimes[it] else it
        infectious <- I(t)
        removed <- R(t)
        plot(object$W, ...)
        title(main = main)
        if (doLegend) do.call(legend, legend.opts)
        if (doTimer) {
            ttxt <- sprintf(timerformat, t)
            do.call(legend, c(list(legend = ttxt), timer.opts))
        }
        if (!is.null(col.influence)) {
            iRids <- which(infectious)
            if (sequential) setdiff(iRids, it)
            for(j in iRids) {
                iR <- shift.owin(object$events@data$.influenceRegion[[j]],
                                 s$eventCoords[j,])
                plot(iR, add = TRUE, col = col.influence, border = NA)
            }
        }
        rTable <- eventCoordsTypes[removed,,drop=FALSE]
        if (nrow(rTable) > 0L) multpoints(rTable, col = col.R)
        iTable <- eventCoordsTypes[infectious,,drop=FALSE]
        if (nrow(iTable) > 0L) multpoints(iTable, col = col.I)
        infectiousNew <- if (sequential) it else infectious & !I(told)
        iTableNew <- eventCoordsTypes[infectiousNew,,drop=FALSE]
        if (nrow(iTableNew) > 0L) multpoints(iTableNew, col = col.current)
        told <- t
        if (verbose) setTxtProgressBar(pb, it)
        Sys.sleep(sleep)
    }
    if (verbose) close(pb)
    invisible(NULL)
}


### plot method for epidataCS

plot.epidataCS <- function (x, aggregate = c("time", "space"), subset, ...)
{
    aggregate <- match.arg(aggregate)
    FUN <- paste("plot.epidataCS", aggregate, sep="_")
    do.call(FUN, args=list(x=quote(x), subset=substitute(subset), ...))
}


### plot.epidataCS(x, aggregate = "time") -> number of cases over time
## in case t0.Date is specified, hist.Date() is used and breaks must set in ... (e.g. "months")

plot.epidataCS_time <- function (x, subset, t0.Date = NULL, freq = TRUE,
    col = "white", add = FALSE,
    xlim = NULL, ylim = NULL, xlab = "Time", ylab = NULL, main = NULL,
    panel.first = abline(h=axTicks(2), lty=2, col="grey"), ...)
{
    timeRange <- with(x$stgrid, c(start[1L], stop[length(stop)]))
    eventTimes <- if (missing(subset)) x$events$time else {
        do.call(base::subset, list(x = quote(marks(x)),
                                   subset = substitute(subset),
                                   select = "time", drop = TRUE))
    }
    if (length(eventTimes) == 0) stop("no events left after 'subset'")
    if (!is.null(t0.Date)) {
        stopifnot(inherits(t0.Date, "Date") || is.vector(t0.Date), length(t0.Date) == 1L)
        t0.Date <- as.Date(t0.Date)
        t0 <- timeRange[1L]
        if (is.null(xlim)) xlim <- t0.Date + (timeRange - t0)
        if (missing(xlab)) xlab <- paste0("Time (", list(...)[["breaks"]], ")")
        eventTimes <- t0.Date + as.integer(eventTimes - t0)
        ## we need integer dates here because otherwise, if the last event
        ## occurs on the last day of a month, year, etc. (depending on
        ## 'breaks') with a fractional date (e.g. as.Date("2009-12-31") + 0.5),
        ## then the automatic 'breaks' (e.g., breaks = "months") will not cover
        ## the data (in the example, it will only reach until
        ## as.Date("2009-12-31")).
    }
    histdata <- if (is.null(t0.Date)) {
        hist(eventTimes, plot=FALSE, warn.unused=FALSE, ...)
    } else {
        hist(eventTimes, plot=FALSE, ...) # warn.unused=FALSE is hard-coded in hist.Date
    }
    if (!add) {
        if (is.null(xlim)) xlim <- timeRange
        if (is.null(ylim)) {
            ylim <- range(0, histdata[[if (freq) "counts" else "density"]])
        }
        if (is.null(ylab)) {
            ylab <- if (freq) "Number of cases" else "Density of cases"
        }
        plot(x=xlim, y=ylim, xlab=xlab, ylab=ylab, main=main, type="n", bty="n")
        force(panel.first)
    }
    plot(histdata, freq = freq, add = TRUE, col = col, ...)
    box()          # because white filling of bars might overdraw the inital box
    invisible(histdata)
}


### plot.epidataCS(x, aggregate = "space") -> spatial point pattern

plot.epidataCS_space <- function (x, subset,
    cex.fun = sqrt, points.args = list(cex=0.5),
    colTypes = rainbow(nlevels(x$events$type)), ...)
{
    stopifnot(is.list(points.args))
    events <- if (missing(subset)) x$events else {
        eval(substitute(base::subset(x$events, subset=.subset),
                        list(.subset=substitute(subset))))
        ## do.call(base::subset,
        ##         args = list(x=quote(x$events), subset=substitute(subset)))
    }
    events@data[["_MULTIPLICITY_"]] <- multiplicity(events)
    events <- events[!duplicated(coordinates(events)),]
    pointcex <- cex.fun(events$"_MULTIPLICITY_")
    pointcex <- pointcex * points.args$cex
    points.args$cex <- NULL
    if (is.null(points.args[["col"]])) {
        points.args$col <- colTypes[events$type]
    }
    plot(x$W, ...)
    do.call("points", c(alist(x=events, cex=pointcex), points.args))
    invisible()
}



######################################################################
# Transform _twinstim_ epidataCS to _twinSIR_ epidata object
######################################################################

# this only generates a SIS epidemic, i.e. atRiskY is set to 1 immediately after recovery
# length of infectious period is taken from epidataCS$events$eps.t
# fcols are not generated here. these must be generated by a second call to twinSIR's as.epidata with desired f. (for safety)
# tileCentroids is a coordinate matrix whose row names are the tile levels
as.epidata.epidataCS <- function (data, tileCentroids, eps = 0.001, ...)
{
    if (!requireNamespace("intervals"))
        stop("conversion from ", dQuote("epidataCS"), " to ", dQuote("epidata"),
             " requires the ", dQuote("intervals"), " package")
    
    ### generate twinSIR's epidata object from stgrid (no events)
    centroidIdx <- match(levels(data$stgrid$tile), rownames(tileCentroids), nomatch = NA_integer_)
    if (any(is.na(centroidIdx))) {
        stop("some levels of 'data$stgrid$tile' are not available from 'tileCentroids'")
    }
    centroids <- tileCentroids[centroidIdx,]
    if (any(c("xCent", "yCent") %in% names(data$stgrid))) {
        stop("'data$stgrid' already has columns \"xCent\" and \"yCent\"")
    }
    stgrid <- cbind(data$stgrid,
        atRiskY = 1L, event = 0L, Revent = 0L,
        xCent = centroids[,1], yCent = centroids[,2]
        # relies on ordering of stgrid by first BLOCK, then tile
    )
    names(stgrid)[names(stgrid)=="tile"] <- "id"

    ### now determine "events" with respect to the tiles
    # individual data
    indItimes <- data$events$time
    if (anyDuplicated(indItimes)) stop("'data$events' has concurrent event times")
    indRtimes <- indItimes + data$events$eps.t
    indInts <- intervals::Intervals(cbind(indItimes, indRtimes, deparse.level = 0L))
    indTiles <- data$events$tile

    # tile data
    tileRows <- tapply(seq_along(indTiles), indTiles, c, simplify = FALSE)
    tileInts <- lapply(tileRows, function (rows) {
        if (length(rows)==0L) { matrix(0,0,2) } else if (length(rows)==1L) {
            as.matrix(indInts[rows])
        } else as.matrix(intervals::reduce(indInts[rows]))
    })
    tileNames <- rep.int(names(tileInts), sapply(tileInts, nrow))
    tileItimes <- unlist(lapply(tileInts, function(ints) ints[,1]), use.names=FALSE)
    tileRtimes <- unlist(lapply(tileInts, function(ints) ints[,2]), use.names=FALSE)

    # there are possibly Rtimes which equal Itimes of other individuals
    # => break ties by considering Rtime shortly before Itime (arbitrary choice)
    while(length(dup <- which(tileRtimes %in% tileItimes)) > 0L) {
        tileRtimes[dup] <- tileRtimes[dup] - eps
    }
    # now there could be duplicated Rtimes... grml (choose another 'eps' in this case)
    if (anyDuplicated(tileRtimes)) {
        stop("breaking ties introduced duplicated Rtimes")
    }

    ### add additional stop times to stgrid for tile infections and recoveries
    requiredStopTimes <- sort(c(tileItimes, tileRtimes))
    class(stgrid) <- c("epidataCS", "data.frame")
    attr(stgrid, "timeRange") <- c(stgrid$start[1], tail(stgrid$stop,1))
    cat("Inserting extra stop times in 'stgrid' (this might take a while)... ")
    evHist <- intersperse(stgrid, requiredStopTimes) # this resets the BLOCK index
    class(evHist) <- "data.frame"
    ### <- THIS IS THE MOST TIME-CONSUMING PART OF THIS FUNCTION !!!
    cat("Done.\n")

    ### set event, Revent and atRiskY indicators
    tileNamesCodes <- match(tileNames, levels(evHist$id))
    # event indicator (currently in evHist event==0 everywhere)
    idxItimes <- match(tileItimes, evHist$stop) - 1L + tileNamesCodes
    evHist$event[idxItimes] <- 1L
    # Revent indicator (currently in evHist Revent==0 everywhere)
    idxRtimes <- match(tileRtimes, evHist$stop) - 1L + tileNamesCodes  # (may contain NA's if Revent after last stop)
    evHist$Revent[idxRtimes] <- 1L
    # atRiskY indicator
    .atRiskY <- rep.int(1L, nrow(evHist))
    nTiles <- nlevels(evHist$id)
    nBlocks <- tail(evHist$BLOCK, 1)
    stopTimes <- unique(evHist$stop)  # has length nBlocks
    for (i in seq_along(tileItimes)) {
        .Itime <- tileItimes[i]
        .Rtime <- tileRtimes[i]
        .tileCode <- tileNamesCodes[i]
        idxsTileInEpi <- seq(.tileCode, by=nTiles, length.out=nBlocks)
        first0block <- match(.Itime, stopTimes) + 1L
        last0block <- if (.Rtime > stopTimes[nBlocks]) nBlocks else match(.Rtime, stopTimes)
        .atRiskY[idxsTileInEpi[first0block:last0block]] <- 0L
    }
    evHist$atRiskY <- .atRiskY

    ### Return final epidata object of twinSIR-type
    cat("Generating final \"epidata\" object for use with twinSIR... ")
    epi <- as.epidata(evHist[-grep("BLOCK", names(evHist))],
        id.col="id", start.col="start", stop.col="stop", atRiskY.col="atRiskY",
        event.col="event", Revent.col="Revent", coords.cols=c("xCent","yCent")
    )
    cat("Done.\n")
    epi
}


####################################################################
### Transform "epidataCS" to "sts" by aggregation of cases on stgrid
####################################################################

epidataCS2sts <- function (object, freq, start,
                           neighbourhood, tiles = NULL,
                           popcol.stgrid = NULL, popdensity = TRUE)
{
    stopifnot(inherits(object, "epidataCS"))
    tileLevels <- levels(object$stgrid$tile)
    if (!is.null(tiles)) {
        stopifnot(inherits(tiles, "SpatialPolygons"),
                  tileLevels %in% row.names(tiles))
        tiles <- tiles[tileLevels,]
    }

    ## prepare sts components
    epoch <- unique(object$stgrid$BLOCK) # epidataCS is sorted
    eventsByCell <- with(object$events@data,
                         table(BLOCK=factor(BLOCK, levels=epoch), tile))
    if (missing(neighbourhood)) { # auto-detect neighbourhood from tiles
        if (is.null(tiles))
            stop("'tiles' is required for auto-generation of 'neighbourhood'")
        neighbourhood <- poly2adjmat(tiles, zero.policy=TRUE)
        if (any(rowSums(neighbourhood) == 0))
            warning("generated neighbourhood matrix contains islands")
    }
    populationFrac <- if (is.null(popcol.stgrid)) NULL else {
        stopifnot(is.vector(popcol.stgrid), length(popcol.stgrid) == 1)
        popByCell <- object$stgrid[[popcol.stgrid]]
        if (popdensity) popByCell <- popByCell * object$stgrid[["area"]]
        totalpop <- sum(popByCell[seq_along(tileLevels)])
        matrix(popByCell/totalpop,
               nrow=length(epoch), ncol=length(tileLevels),
               byrow=TRUE, dimnames=dimnames(eventsByCell))
    }

    ## initialize sts object
    new("sts", epoch = epoch, freq=freq, start=start,
        observed=unclass(eventsByCell), neighbourhood=neighbourhood,
        populationFrac=populationFrac, map=tiles, epochAsDate=FALSE)
}


###################################################
### Distances from potential (eps.s, eps.t) sources
###################################################

getSourceDists <- function (object)
{
    ## extract required info from "epidataCS"-object
    eventCoords <- coordinates(object$events)
    .sources <- object$events$.sources

    ## distance matrix and number of sources
    distmat <- as.matrix(dist(eventCoords))
    nsources <- sapply(.sources, length)
    hasSources <- nsources > 0
    cnsources <- c(0, cumsum(nsources))

    ## generate vector of distances of events to their potential sources
    sourcedists <- numeric(sum(nsources))
    for (i in which(hasSources)) {
        .sourcesi <- .sources[[i]]
        .sourcedists <- distmat[i, .sourcesi]
        .idx <- cnsources[i] + seq_len(nsources[i])
        sourcedists[.idx] <- .sourcedists
        names(sourcedists)[.idx] <- paste(i, .sourcesi, sep="<-")
    }

    ## Done
    sourcedists
}
