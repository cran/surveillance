################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Formulation of an endemic-only twinstim as a pseudo Poisson-GLM
###
### Copyright (C) 2013 Sebastian Meyer
### $Revision: 484 $
### $Date: 2013-01-24 16:51:02 +0100 (Do, 24. Jan 2013) $
################################################################################


### a type-independent endemic-only twinstim equals a weighted Poisson-GLM on
### the stgrid data with outcome nEvents/weights, and weights = K * dt * ds

glm.epidataCS <- function (data, formula = NULL)
{
    if (is.null(formula)) {
        formula <- as.formula(paste("~",
                       paste(names(data$stgrid)[-c(1,3:5)], collapse="+")
                   ))
    } else if ("1 | type" %in% attr(terms(formula), "term.labels")) {
        stop("type-specific endemic intercepts are currently not supported")
        ## if we had a type-specific model, we would have to set up a
        ## "stkappagrid", i.e. with nBlocks*nTiles*nTypes rows, and remove the
        ## factor nTypes from the weights
    }
    
    ## aggregated number of events in each cell of the stgrid
    eventsbycell <- c(table(with(data$events@data,
                                 interaction(tile, BLOCK,
                                             drop=FALSE, sep=".", lex.order=FALSE)
                                 )))
    nEvents <- eventsbycell[paste(data$stgrid$tile, data$stgrid$BLOCK, sep=".")]
    nEvents[is.na(nEvents)] <- 0L
    stopifnot(sum(nEvents) == nrow(data$events))

    ## weights
    nTypes <- nlevels(data$events$type)
    weights <- with(data$stgrid, nTypes * (stop-start) * area)

    ## Poisson-GLM
    environment(formula) <- environment() # -> nEvents and weights are visible
    glm(update(formula, nEvents/weights ~ .), data=data$stgrid,
        family=poisson, weights=weights)
}
