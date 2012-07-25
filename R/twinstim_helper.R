################################################################################
### Some internal helper functions for 'twinstim'.
###
### Author: Sebastian Meyer
### $Date: 2010-11-16 04:08:22 +0100 (Tue, 16 Nov 2010) $
################################################################################


## Compute distance from points to boundary
## (adapted from spatstat::bdist.points, DEPENDS ON spatstat::distppl)

# xy is the coordinate matrix of the points
# poly is a polygonal domain of class "gpc.poly"
# the function does not check if points are actually inside the polygonal domain
bdist <- function (xy, poly)
{
    result <- rep(Inf, nrow(xy))
    bdry <- poly@pts
    for (i in seq_along(bdry)) {
        polly <- bdry[[i]]
        nsegs <- length(polly$x)
        for (j in 1L:nsegs) {
            j1 <- if (j < nsegs) j + 1L else 1L
            seg <- c(polly$x[j], polly$y[j], polly$x[j1], polly$y[j1])
            result <- pmin(result, spatstat::distppl(xy, seg))
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



