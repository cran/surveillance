################################################################################
### Wrapper function for fanplot::fan()
###
### Copyright (C) 2017-2018 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

fanplot <- function (quantiles, probs, means = NULL, observed = NULL,
    start = 1, fan.args = list(), means.args = list(), observed.args = list(),
    key.args = NULL, xlim = NULL, ylim = NULL,
    xlab = "Time", ylab = "No. infected", add = FALSE, ...)
{
    stopifnot(is.matrix(quantiles), length(probs) == ncol(quantiles),
              is.null(means) || length(means) == nrow(quantiles),
              is.null(observed) || length(observed) == nrow(quantiles),
              isScalar(start))

    ## axis range
    if (is.null(xlim))
        xlim <- c(1 - 0.5, nrow(quantiles) + 0.5) + (start-1)
    if (is.null(ylim))
        ylim <- c(0, max(quantiles, observed))

    ## graphical parameters
    stopifnot(is.list(fan.args))
    fan.args <- modifyList(
        list(data = t(quantiles), data.type = "values", probs = probs,
             start = start, fan.col = heat.colors, ln = NULL),
        fan.args, keep.null = TRUE)

    ## initialize empty plot
    if (!add)
        plot(xlim, ylim, type = "n", xlab = xlab, ylab = ylab, ...)

    ## add fan
    do.call(fanplot::fan, fan.args)

    ## add point predictions
    if (!is.null(means) && is.list(means.args)) {
        means.args <- modifyList(
            list(x = seq_along(means) + (start-1), y = means,
                 type = "l", lwd = 2, col = "white"),
            means.args)
        do.call("lines", means.args)
    }

    ## add observed time series
    if (!is.null(observed) && is.list(observed.args)) {
        observed.args <- modifyList(
            list(x = seq_along(observed) + (start-1), y = observed,
                 type = "b", lwd = 2),
            observed.args)
        do.call("lines", observed.args)
    }

    ## add color key
    if (is.list(key.args)) {
        key.args <- modifyList(
            list(start = xlim[2L] - 1, ylim = c(ylim[1L] + mean(ylim), ylim[2L]),
                 data.type = "values", style = "boxfan", probs = fan.args$probs,
                 fan.col = fan.args$fan.col, ln = NULL, space = 0.9,
                 rlab = quantile(fan.args$probs, names = FALSE, type = 1)),
            key.args)
        ## convert ylim to data
        key.args$data <- matrix(seq.int(from = key.args$ylim[1L], to = key.args$ylim[2L],
                                        length.out = length(fan.args$probs)))
        key.args$ylim <- NULL
        tryCatch(do.call(fanplot::fan, key.args), error = function (e)
            warning("color key could not be drawn, probably due to non-standard 'probs'",
                    call. = FALSE))
    }

    invisible(NULL)
}
