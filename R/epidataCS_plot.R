################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### plot-method for "epidataCS" objects
###
### Copyright (C) 2009-2014 Sebastian Meyer
### $Revision: 903 $
### $Date: 2014-04-10 22:18:03 +0200 (Thu, 10 Apr 2014) $
################################################################################


plot.epidataCS <- function (x, aggregate = c("time", "space"), subset, ...)
{
    aggregate <- match.arg(aggregate)
    FUN <- paste("epidataCSplot", aggregate, sep="_")
    do.call(FUN, args=list(x=quote(x), subset=substitute(subset), ...))
}


### plot.epidataCS(x, aggregate = "time") -> number of cases over time
## in case t0.Date is specified, hist.Date() is used and breaks must set in ... (e.g. "months")

epidataCSplot_time <- function (x, subset, t0.Date = NULL, freq = TRUE,
    col = rainbow(nTypes), cumulative = list(), add = FALSE, mar = NULL,
    xlim = NULL, ylim = NULL, xlab = "Time", ylab = NULL, main = NULL,
    panel.first = abline(h=axTicks(2), lty=2, col="grey"),
    legend.types = list(), ...)
{
    timeRange <- with(x$stgrid, c(start[1L], stop[length(stop)]))
    eventTimesTypes <- if (missing(subset)) {
        x$events@data[c("time", "type")]
    } else {
        do.call(base::subset, list(x = quote(marks.epidataCS(x)),
                                   subset = substitute(subset),
                                   select = c("time", "type")))
    }
    if (nrow(eventTimesTypes) == 0L) stop("no events left after 'subset'")
    if (is.list(cumulative)) {
        csums <- tapply(eventTimesTypes$time, eventTimesTypes["type"],
                        function (t) cumsum(table(t)), simplify=FALSE)
        if (is.null(cumulative[["axis"]])) cumulative[["axis"]] <- TRUE
    }
    typeNames <- levels(eventTimesTypes$type)
    nTypes <- length(typeNames)
    eventTimesTypes$type <- as.integer(eventTimesTypes$type)
    col <- rep_len(col, nTypes)
    
    if (!is.null(t0.Date)) {
        stopifnot(length(t0.Date) == 1L)
        t0.Date <- as.Date(t0.Date)
        t0 <- timeRange[1L]
        if (is.null(xlim)) xlim <- t0.Date + (timeRange - t0)
        if (missing(xlab)) xlab <- paste0("Time (", list(...)[["breaks"]], ")")
        eventTimesTypes$time <- t0.Date + as.integer(eventTimesTypes$time - t0)
        ## we need integer dates here because otherwise, if the last event
        ## occurs on the last day of a month, year, etc. (depending on
        ## 'breaks') with a fractional date (e.g. as.Date("2009-12-31") + 0.5),
        ## then the automatic 'breaks' (e.g., breaks = "months") will not cover
        ## the data (in the example, it will only reach until
        ## as.Date("2009-12-31")).
    }
    gethistdata <- function (types = seq_len(nTypes)) {
        times <- eventTimesTypes$time[eventTimesTypes$type %in% types]
        if (is.null(t0.Date)) {
            hist(times, plot=FALSE, warn.unused=FALSE, ...)
        } else {
            hist(times, plot=FALSE, ...)
            ## warn.unused=FALSE is hard-coded in hist.Date
        }
    }
    histdata <- gethistdata()
    if (!add) {
        if (is.null(xlim)) xlim <- timeRange
        if (is.null(ylim)) {
            ylim <- range(0, histdata[[if (freq) "counts" else "density"]])
        }
        if (is.null(ylab)) {
            ylab <- if (freq) "Number of cases" else "Density of cases"
        }
        
        if (is.null(mar)) {
            mar <- par("mar")
            if (is.list(cumulative) && cumulative$axis) mar[4L] <- mar[2L]
        }
        opar <- par(mar = mar); on.exit(par(opar))
        plot(x=xlim, y=ylim, xlab=xlab, ylab=ylab, main=main, type="n", bty="n")
        force(panel.first)
    }

    ## plot histogram (over all types)
    plot(histdata, freq = freq, add = TRUE, col = col[1L], ...)
    box()          # because white filling of bars might overdraw the inital box

    ## add type-specific sub-histograms
    typesEffective <- sort(unique(eventTimesTypes$type))
    for (typeIdx in seq_along(typesEffective)[-1L]) {
        plot(gethistdata(typesEffective[typeIdx:length(typesEffective)]),
             freq = freq, add = TRUE, col = col[typesEffective[typeIdx]], ...)
    }

    ## optionally add cumulative number of cases
    if (is.list(cumulative)) {
        aT2 <- axTicks(2)
        div <- length(aT2) - 1L
        cumulative <- modifyList(
            list(maxat = ceiling(max(unlist(csums))/div)*div,
                 col = apply(col2rgb(col)/255*0.6, 2, function (x)
                             rgb(x[1L], x[2L], x[3L])),
                 lwd = 3, axis = TRUE, lab = "Cumulative number of cases"),
            cumulative)
        csum2y <- function (x) x / cumulative$maxat * aT2[length(aT2)]
        for (typeIdx in typesEffective) {
            .times <- as.numeric(names(csums[[typeIdx]]))
            lines(if (is.null(t0.Date)) .times else t0.Date + .times - t0,
                  csum2y(csums[[typeIdx]]), lwd=cumulative$lwd,
                  col=cumulative$col[typeIdx])
        }
        if (cumulative$axis) {
            axis(4, at=aT2, labels=aT2/aT2[length(aT2)]*cumulative$maxat)
            mtext(cumulative$lab, side=4, line=3, las=0)
        }
    }

    ## optionally add legend
    if (is.list(legend.types) && length(typesEffective) > 1) {
        legend.types <- modifyList(
            list(x="topleft", legend=typeNames[typesEffective],
                 title="Type", fill=col[typesEffective]),
            legend.types)
        do.call("legend", legend.types)
    }
    
    invisible(histdata)
}


### plot.epidataCS(x, aggregate = "space") -> spatial point pattern

epidataCSplot_space <- function (x, subset,
    cex.fun = sqrt, points.args = list(), add = FALSE,
    legend.types = list(), legend.counts = list(), ...)
{
    events <- if (missing(subset)) x$events else {
        eval(substitute(base::subset(x$events, subset=.subset),
                        list(.subset=substitute(subset))))
    }
    eventCoordsTypesCounts <- countunique(
        cbind(coordinates(events), type = as.integer(events$type))
        )
    pointCounts <- eventCoordsTypesCounts[,"COUNT"]
    typeNames <- levels(events$type)
    nTypes <- length(typeNames)
    typesEffective <- sort(unique(eventCoordsTypesCounts[,"type"]))

    ## point style
    colTypes <- list(...)[["colTypes"]]  # backwards compatibility for < 1.8
    if (is.null(colTypes)) {
        colTypes <- rainbow(nTypes)
    } else warning("argument 'colTypes' is deprecated; ",
                   "use 'points.args$col' instead")
    points.args <- modifyList(list(pch=1, col=colTypes, lwd=1, cex=0.5),
                              points.args)
    styleArgs <- c("pch", "col", "lwd")
    points.args[styleArgs] <- lapply(points.args[styleArgs],
                                     rep_len, length.out=nTypes)

    ## select style parameters according to the events' types
    points.args_pointwise <- points.args
    points.args_pointwise[styleArgs] <- lapply(
        points.args_pointwise[styleArgs], "[",
        eventCoordsTypesCounts[,"type"])
    points.args_pointwise$cex <- points.args_pointwise$cex * cex.fun(pointCounts)
    
    ## plot
    if (!add) plot(x$W, ...)
    do.call("points", c(alist(x=eventCoordsTypesCounts[,1:2]),
                        points.args_pointwise))

    ## optionally add legends
    if (is.list(legend.types) && length(typesEffective) > 1) {
        legend.types <- modifyList(
            list(x="topright", legend=typeNames[typesEffective],
                 title="Type",
                 #pt.cex=points.args$cex, # better use par("cex")
                 pch=points.args$pch[typesEffective],
                 col=points.args$col[typesEffective],
                 pt.lwd=points.args$lwd[typesEffective]),
            legend.types)
        do.call("legend", legend.types)
    }
    if (is.list(legend.counts) && any(pointCounts > 1)) {
        if (is.null(legend.counts[["counts"]])) {
            counts <- unique(round(10^(do.call("seq", c(
                as.list(log10(range(pointCounts))), list(length.out=5)
                )))))
        } else {
            counts <- as.vector(legend.counts[["counts"]], mode="integer")
            legend.counts[["counts"]] <- NULL
        }
        legend.counts <- modifyList(
            list(x="bottomright", bty="n", legend=counts,
                 pt.cex=points.args$cex * cex.fun(counts),
                 pch=points.args$pch[1L],
                 col=if(length(unique(points.args$col)) == 1L)
                 points.args$col[1L] else 1,
                 pt.lwd=points.args$lwd[1L]),
            legend.counts)
        do.call("legend", legend.counts)
    }
    invisible()
}
