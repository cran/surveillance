################################################################################
### calibrationTest() for "hhh4" fits
###
### Copyright (C) 2015,2017 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################

calibrationTest.hhh4 <- function (x,
                                  subset = x$control$subset,
                                  units = seq_len(x$nUnit),
                                  ...)
{
    ## perform the calibration test in the specified subset
    res <- calibrationTest.default(
        x = x$stsObj@observed[subset, units, drop = FALSE],
        mu = x$fitted.values[match(subset, x$control$subset), units, drop = FALSE],
        size = psi2size.hhh4(x, subset, units),
        ...)

    ## change "data.name" to be the name of the supplied model
    res$data.name <- deparse(substitute(x))
    res
}

calibrationTest.oneStepAhead <- function (x, units = NULL, ...)
{
    ## perform the calibration test
    res <- if (is.null(units)) {
        calibrationTest.default(
            x = x$observed,
            mu = x$pred,
            size = psi2size.oneStepAhead(x),
            ...)
    } else {
        calibrationTest.default(
            x = x$observed[, units, drop = FALSE],
            mu = x$pred[, units, drop = FALSE],
            size = psi2size.oneStepAhead(x)[, units, drop = FALSE],
            ...)
    }

    ## change "data.name" to be the name of the supplied "oneStepAhead" object
    res$data.name <- deparse(substitute(x))
    res

}
