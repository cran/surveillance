################################################################################
### Plot a surveillance time series ("sts") object using ggplot2
###
### Copyright (C) 2018 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

autoplot.sts <- function (object, population = FALSE,
                          units = NULL, as.one = FALSE,
                          scales = "fixed", ...)
{
    stopifnot(is(object, "sts"))
    data <- tidy.sts(object)

    ## select subset of units to plot
    if (!is.null(units)) {
        ## ensure that 'units' are labels, not indices
        units <- unname(setNames(nm = levels(data$unit))[units])
        data <- data[data$unit %in% units, , drop=FALSE]
    }

    ## scale counts by population
    if (doInc <- isScalar(population) || isTRUE(population))
        data$observed <- data$observed / (data$population / population)

    p <- ggplot2::ggplot(
        data = data,
        mapping = ggplot2::aes_(x = ~date, y = ~observed, group = ~unit)
    )
    if (as.one) {
        p <- p + ggplot2::geom_line(ggplot2::aes_(colour = ~unit))
    } else {
        p <- p + ggplot2::geom_bar(stat = "identity") +
            ggplot2::facet_wrap(~unit, scales = scales, drop = TRUE)
    }
    p + ggplot2::labs(x = "Time", y = if(doInc) "Incidence" else "No. infected")
}
