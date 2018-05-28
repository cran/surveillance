################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Animated map (and time series chart) of an sts-object (or matrix of counts)
###
### Copyright (C) 2013-2016,2018 Sebastian Meyer
### $Revision: 2123 $
### $Date: 2018-05-03 16:29:46 +0200 (Thu, 03. May 2018) $
################################################################################


### Corresponding to the S3-generic function animate(),
### we define a method for the S4-class "sts" and omit the recommended
### setGeneric("animate"); setMethod("animate", "sts", animate.sts)
### [see Section "Methods for S3 Generic Functions" in help("Methods")]

animate.sts <- function (object, tps = NULL, cumulative = FALSE,
                         population = NULL, at = 10, ...,
                         timeplot = list(height = 0.3, fill = FALSE),
                         sleep = 0.5, verbose = interactive(), draw = TRUE)
{
    if (draw && dev.interactive())
        message("Advice: use facilities of the \"animation\" package, e.g.,\n",
                "        saveHTML() to view the animation in a web browser.")

    if (is.null(tps))
        tps <- seq_len(nrow(object))
    if (!is.null(population)) { # get population matrix
        population <- parse_population_argument(population, object)
    }

    ## determine color breaks (checkat() is defined in stsplot_space.R)
    at <- checkat(at, data=.rangeOfDataToPlot(object, tps, cumulative, population),
                  counts=is.null(population))

    ## style of the additional temporal plot
    if (is.list(timeplot)) {
        timeplot <- modifyList(eval(formals()$timeplot), timeplot)
        timeplot_height <- timeplot$height
        timeplot_fill <- timeplot$fill
        timeplot$height <- timeplot$fill <- NULL  # not for stsplot_timeSimple()
        stopifnot(timeplot_height > 0, timeplot_height < 1)
    }

    if (verbose)
        pb <- txtProgressBar(min=0, max=length(tps), initial=0, style=3)
    grobs <- vector(mode = "list", length = length(tps))
    for(i in seq_along(tps)) {
        cti <- if (cumulative) seq_len(i) else i
        ls <- stsplot_space(object, tps=tps[cti], population=population,
                            at=at, ...)
        if (is.list(timeplot) && requireNamespace("gridExtra")) {
            stopifnot(packageVersion("gridExtra") >= "2.0.0")
            lt <- do.call("stsplot_timeSimple", c(
                list(x=object, tps=tps, highlight=cti),
                timeplot))
            if (!timeplot_fill) {
                lt$aspect.fill <- FALSE
                lt$aspect.ratio <- timeplot_height * ls$aspect.ratio
            }
            grobs[[i]] <- gridExtra::arrangeGrob(
                ls, lt, heights=c(1-timeplot_height, timeplot_height))
            ## alternative using package "gtable":
            ## drawDetails.lattice <- function (x, recording = FALSE)
            ##     plot(x$p, newpage = FALSE)
            ## heights <- c(1-timeplot_height, timeplot_height)
            ## gt <- gtable::gtable(widths = grid::unit(1, units = "null"),
            ##                      heights = grid::unit(heights, units = "null"))
            ## gt <- gtable::gtable_add_grob(gt, list(grid::grob(p = ls, cl = "lattice"),
            ##                                        grid::grob(p = lt, cl = "lattice")),
            ##                               t = 1:2, l = 1)
            if (draw) {
                grid::grid.newpage()
                grid::grid.draw(grobs[[i]])
            }
        } else {
            grobs[[i]] <- ls
            if (draw) print(ls)
        }
        if (verbose) setTxtProgressBar(pb, i)
        if (dev.interactive()) Sys.sleep(sleep)
    }
    if (verbose) close(pb)
    invisible(grobs)
}


### additional time plot below the map

stsplot_timeSimple <- function (x, tps = NULL, highlight = integer(0),
                                inactive = list(col="gray", lwd=1),
                                active = list(col=1, lwd=4),
                                as.Date = x@epochAsDate, ...)
{
    observed <- if (inherits(x, "sts")) observed(x) else x
    if (is.null(tps)) {
        tps <- seq_len(nrow(observed))
    } else {
        observed <- observed[tps,,drop=FALSE]
    }
    epoch <- if (inherits(x, "sts")) epoch(x, as.Date = as.Date)[tps] else tps

    if (anyNA(observed))
        warning("ignoring NA counts in time series plot")

    ## build highlight-specific style vectors (col, lwd, ...)
    stopifnot(is.list(inactive), is.list(active))
    stylepars <- intersect(names(inactive), names(active))
    styleargs <- sapply(stylepars, function (argname) {
        res <- rep.int(inactive[[argname]], length(tps))
        res[highlight] <- active[[argname]]
        res
    }, simplify=FALSE, USE.NAMES=TRUE)

    par_no_top_padding <- list(
        layout.heights = list(top.padding = 0,
                              main.key.padding = 0,
                              key.axis.padding = 0)
    )
    xyplot.args <- modifyList(
        c(list(x = rowSums(observed, na.rm = TRUE) ~ epoch,
               type = "h", ylab = "", xlab = "",
               par.settings = par_no_top_padding),
          styleargs),
        list(...))
    do.call(lattice::xyplot, xyplot.args)
}


### determine data range for automatic color breaks 'at'

.rangeOfDataToPlot <- function (object, tps, cumulative = FALSE,
                                population = NULL)
{
    observed <- if (inherits(object, "sts")) observed(object) else object
    observed <- observed[tps,,drop=FALSE]
    if (!is.null(population)) { # compute (cumulative) incidence
        observed <- if (cumulative) {
            observed / rep(population[tps[1L],], each = nrow(observed))
        } else {
            observed / population[tps,,drop=FALSE]
        }
    }
    range(if (cumulative) c(observed[1L,], colSums(observed)) else observed,
          na.rm = TRUE)
}
