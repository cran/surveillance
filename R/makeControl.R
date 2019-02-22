################################################################################
### Convenient construction of a list of control arguments for "hhh4" models
###
### Copyright (C) 2014-2015 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

makeControl <- function (f = list(~1), S = list(0, 0, 1), period = 52,
                         offset = 1, ...)
{
    ## set model components
    control <- mapply(function (f, S, period, offset) {
        f <- addSeason2formula(f = f, S = S, period = period)
        list(f = f, offset = offset)
    }, f, S, period, offset, SIMPLIFY = FALSE, USE.NAMES = FALSE)
    names(control) <- c("ar", "ne", "end")
    ## default: negative-binomial distribution with common overdispersion
    control$family <- "NegBin1"
    ## customization via ... arguments
    modifyList(control, list(...))
}
