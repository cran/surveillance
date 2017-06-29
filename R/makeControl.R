################################################################################
### Convenient construction of a list of control arguments for "hhh4" models
###
### Copyright (C) 2014-2015 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
################################################################################

##' Generate \code{control} Settings for an \code{hhh4} Model
##'
##' @param f,S,period arguments for \code{\link{addSeason2formula}} defining
##'     each of the three model formulae in the order (\code{ar}, \code{ne},
##'     \code{end}). Recycled if necessary within \code{\link{mapply}}.
##' @param offset multiplicative component offsets in the order (\code{ar},
##'     \code{ne}, \code{end}).
##' @param ... further elements for the \code{\link{hhh4}} control list. The
##'     \code{family} parameter is set to \code{"NegBin1"} by default.
##' @return a list for use as the \code{control} argument in \code{\link{hhh4}}.
##' @author Sebastian Meyer
##' @examples
##' makeControl()
##'
##' ## a simplistic model for the fluBYBW data
##' ## (first-order transmission only, no district-specific intercepts)
##' data("fluBYBW")
##' mycontrol <- makeControl(
##'     f = list(~1, ~1, ~t), S = c(1, 1, 3),
##'     offset = list(population(fluBYBW)),  # recycled -> in all components
##'     ne = list(normalize = TRUE), verbose = TRUE)
##' str(mycontrol)
##' \dontrun{fit <- hhh4(fluBYBW, mycontrol)}
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
