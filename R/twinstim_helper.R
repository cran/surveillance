################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Some internal helper functions for "twinstim".
###
### Copyright (C) 2009-2012 Sebastian Meyer
### $Revision: 446 $
### $Date: 2012-10-23 15:54:41 +0200 (Di, 23. Okt 2012) $
################################################################################


## Compute distance from points to boundary
## copied in part from function bdist.points() of the "spatstat" package authored
## by A. Baddeley and R. Turner (DEPENDS ON spatstat::distppl)

# xy is the coordinate matrix of the points
# poly is a polygonal domain of class "gpc.poly"
# the function does not check if points are actually inside the polygonal domain
bdist <- function (xy, poly)
{
    result <- rep(Inf, nrow(xy))
    bdry <- poly@pts
    for (i in seq_along(bdry)) {
        polly <- bdry[[i]]
        px <- polly$x
        py <- polly$y
        nsegs <- length(px)
        for (j in seq_len(nsegs)) {
            j1 <- if (j < nsegs) j + 1L else 1L
            seg <- c(px[j], py[j], px[j1], py[j1])
            result <- pmin(result, distppl(xy, seg))
        }
    }
    return(result)
}



### Determines indexes of potential sources of infection

## determine potential sources of the i'th event
## all arguments but i and qmatrix are nEvents-vectors
## -> determine potential sources for eventTimes[i], eventsTypes[i] with
## distances distvec_j = ||s_i - s_j||
determineSources <- function (i, eventTimes, removalTimes, distvec, eps.s,
    eventTypes = NULL, qmatrix)
{
    tp <- eventTimes[i]
    type <- eventTypes[i]   # NULL[i] -> NULL
    infectivity <- (eventTimes < tp) & (removalTimes >= tp)
    #<- eventTimes<tp, not "=" because CIF is left-continuous.
    #Also guarantees no self-infection
    proximity <- as.vector(distvec <= eps.s, mode = "logical")
    #<- as.vector to remove names
    matchType <- if (is.null(eventTypes)) TRUE else {
        typeInfective <- qmatrix[,type] # indexed by integer code of factor
        #<- logical vector indicating for each type if it could infect type of i
        as.vector(typeInfective, mode = "logical")[eventTypes]
        #<- as.vector to remove names
    }
    # return IDs of potential epidemic sources
    which(infectivity & proximity & matchType)
}

## determine the .sources for an epidataCS object, i.e.
## lapply the previous function to all of object$events
determineSources.epidataCS <- function (object)
{
    eventTimes <- object$events$time
    removalTimes <- eventTimes + object$events$eps.t
    eventDists <- as.matrix(dist(object$events@coords, method = "euclidean"))
    lapply(seq_along(eventTimes), function (i) {
        determineSources(i, eventTimes, removalTimes, eventDists[i,],
                         object$events$eps.s, object$events$type,
                         object$qmatrix) 
    })
}



### Check matrix Q

checkQ <- function (qmatrix, typeNames)
{
    nTypes <- length(typeNames)
    qmatrix <- as.matrix(qmatrix)
    stopifnot(nrow(qmatrix) == ncol(qmatrix))
    if (is.null(dimnames(qmatrix))) {
        if (nrow(qmatrix) != nTypes) {
            stop("'qmatrix' must be a ", nTypes, "x", nTypes, " matrix")
        }
        dimnames(qmatrix) <- list(typeNames, typeNames)
    } else {
        stopifnot(rownames(qmatrix) == colnames(qmatrix))
        typesIdx <- match(typeNames, rownames(qmatrix), nomatch=NA_integer_)
        if (idx <- match(TRUE, is.na(typesIdx), nomatch=0L)) {
            stop("missing entries for type '", typeNames[idx], "' in 'qmatrix'")
        }
        qmatrix <- qmatrix[typesIdx,typesIdx,drop=FALSE]
    }
    storage.mode(qmatrix) <- "logical"   # convert entries to TRUE/FALSE by as.logical
    qmatrix
}


### Get row index of 'stgrid' where an event is located (spatio-temporally)
### Here, search BLOCK such that t in (start;stop], i.e. an event at 'stop' is
### still attributed to the previous interval

gridcellOfEvent <- function (t, tilename, stgrid)
{
    ## idx <- with(stgrid, which(tile == tilename & start < t & stop >= t))
    
    ## ~5x faster alternative assuming a full BLOCK x tile grid, which is
    ## sorted by BLOCK and tile (tile varying first), specifically there must be
    ## all levels(stgrid$tile) in every BLOCK in that order;
    ## this structure is guaranteed by checkstgrid()
    blockstart <- match(TRUE, stgrid$stop >= t)
    idx <- blockstart + match(tilename, levels(stgrid$tile)) - 1L
    
    lidx <- length(idx)
    if (lidx == 0L) NA_integer_ else if (lidx == 1L) idx else {
        stop("'stgrid' has overlapping spatio-temporal grid cells")
    }
}



### Really internal helper function, which constructs the function that
### integrates the two-dimensional 'siaf' function over the influence regions of
### the events. The only argument of the returned function is 'siafpars'.
### The returned function is defined in the callers environment, where the
### variables used in the function are available (inside twinstim() or
### simEpidataCS()).

.siafIntFUN <- function (siaf,
   noCircularIR #= all(eps.s > bdist) = all(sapply(influenceRegion, function(x)
                #                           is.null(attr(x,"radius"))))
){
    ## the following variables are unused here, because the environment of
    ## FUN will be set to the parent.frame(), where the variables exist
    ## they are only included to avoid the notes in R CMD check 
    iRareas <- influenceRegion <- eventTypes <- eps.s <- bdist <- nTypes <- NULL
    
    ## define the siaf integration function depending on the siaf specification 
    FUN <- if (attr(siaf, "constant"))
    {
        if (exists("iRareas", where=parent.frame(), mode="numeric")) {
            ## in twinstim(), 'iRareas' are pre-defined to save
            ## computation time (data are fixed during fitting)
            function (siafpars) iRareas
        } else {
            function (siafpars) sapply(influenceRegion, attr, "area")
        }
    } else if (is.null(siaf$Fcircle) || # if siaf$Fcircle not available
               (is.null(siaf$effRange) && noCircularIR))
    {
        function (siafpars, ...) {
            siafInts <- sapply(seq_along(influenceRegion), function (i)
                siaf$F(influenceRegion[[i]], siaf$f, siafpars, eventTypes[i], ...)
            )
            siafInts
        }
    } else if (is.null(siaf$effRange)) # use Fcircle but only delta-trick
    {
        function (siafpars, ...) {
            ## Compute the integral of 'siaf' over each influence region
            siafInts <- numeric(length(influenceRegion))
            for(i in seq_along(siafInts)) {
                eps <- eps.s[i]
                bdisti <- bdist[i]
                siafInts[i] <- if (eps <= bdisti) {
                    ## influence region is completely inside W
                    siaf$Fcircle(eps, siafpars, eventTypes[i])
                } else {
                    ## numerically integrate over polygonal influence region
                    siaf$F(influenceRegion[[i]], siaf$f, siafpars,
                           eventTypes[i], ...)
                }
            }
            siafInts
        }
    } else { # fast Fcircle integration considering the delta-trick AND effRange
        .ret <- function (siafpars, ...) {
            ## Compute computationally effective range of the 'siaf' function
            ## for the current 'siafpars' for each event (type)
            effRangeTypes <- rep(siaf$effRange(siafpars), length.out=nTypes)
            effRanges <- effRangeTypes[eventTypes]   # N-vector
            ## Compute the integral of 'siaf' over each influence region
            siafInts <- numeric(length(influenceRegion))
            for(i in seq_along(siafInts)) {
                eps <- eps.s[i]
                bdisti <- bdist[i]
                effRange <- effRanges[i]
                siafInts[i] <- if (eps <= bdisti) { # influence region is completely inside W
                    siaf$Fcircle(eps, siafpars, eventTypes[i])
                } else if (effRange <= bdisti) { # effective region completely inside W
                    siaf$Fcircle(bdisti, siafpars, eventTypes[i])
                } else { # integrate over polygonal influence region
                    siaf$F(influenceRegion[[i]], siaf$f, siafpars,
                           eventTypes[i], ...)
                }
            }
            siafInts
        }
        if (exists("effRangeTypes", where=parent.frame(), mode="numeric")) {
            ## in simEpidataCS effRangeTypes is pre-calculated outside siafInt to
            ## save computation time ('siafpars' is constant during simulation)
            body(.ret)[[grep("^effRangeTypes <-", body(.ret))]] <- NULL
        }
        .ret
    }
    
    ## set the environment of the siafInt function to the callers environment
    ## (i.e. inside twinstim() or simEpidataCS())
    ## where the variables used in the function are defined
    environment(FUN) <- parent.frame()
    FUN
}


### Helper function, which constructs the function that integrates the 'tiaf'.
### The returned function is defined in the callers environment, where the
### variables used in the function are available (inside twinstim() or
### simEpidataCS()).

.tiafIntFUN <- function ()
{
    ## the following variables are unused here, because the environment of
    ## FUN will be set to the parent.frame(), where the variables exist
    ## they are only included to avoid the notes in R CMD check 
    gIntLower <- gIntUpper <- eventTypes <- tiaf <- NULL
    
    ## from, to and type may be vectors of compatible lengths
    FUN <- function(tiafpars, from = gIntLower, to = gIntUpper,
                    type = eventTypes, G = tiaf$G)
    {
        tiafIntUpper <- G(to, tiafpars, type)
        tiafIntLower <- G(from, tiafpars, type)
        tiafIntUpper - tiafIntLower
    }
    
    ## set the environment of the tiafInt function to the callers environment
    ## (i.e. inside twinstim() or simEpidataCS())
    ## where the default argument values are defined
    environment(FUN) <- parent.frame()
    FUN
}


### rename control arguments with optim names to have names compatible with nlminb

control2nlminb <- function (control, defaults)
{
    renamelist <- cbind(optim  = c("maxit", "REPORT", "abstol", "reltol"),
                        nlminb = c("iter.max", "trace", "abs.tol", "rel.tol"))
    for (i in which(renamelist[,"optim"] %in% names(control))) {
        fromname <- renamelist[i, "optim"]
        toname <- renamelist[i, "nlminb"]
        if (is.null(control[[toname]])) {
            control[[toname]] <- control[[fromname]]
        }
        control[[fromname]] <- NULL
    }
    defaults[names(control)] <- control
    defaults
}
