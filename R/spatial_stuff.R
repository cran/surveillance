################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Spatial helper functions
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 631 $
### $Date: 2013-08-28 17:52:52 +0200 (Mit, 28 Aug 2013) $
################################################################################


### Compute distance from points to boundary
### copied in part from function bdist.points() of the "spatstat" package authored
### by A. Baddeley and R. Turner (DEPENDS ON spatstat::distppl)

## xy is the coordinate matrix of the points
## poly is a polygonal domain of class "gpc.poly"
## the function does not check if points are actually inside the polygonal domain
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
            ## calculate distances of points to segment j of polygon i
            dists <- distppl(xy, seg)
            result <- pmin(result, dists)
        }
    }
    return(result)
}


### sample n points uniformly on a disc with radius r

runifdisc <- function (n, r = 1)
{
    rangle <- runif(n, 0, 2*pi)
    rdist <- r * sqrt(runif(n, 0, 1))
    rdist * cbind(cos(rangle), sin(rangle))
}


### Redefinition of gpclib's scale.poly method to also do centering

scale.gpc.poly <-
    function (x, center = c(0,0), scale = c(1,1)) {
        x@pts <- lapply(x@pts, function (p) {
            p$x <- (p$x-center[1]) / scale[1]
            p$y <- (p$y-center[2]) / scale[2]
            p
        })
        x
    }


### Same as inside.owin for gpc.poly (using point.in.polygon from package sp)

inside.gpc.poly <- function(x, y = NULL, polyregion, mode.checked = FALSE)
{
    xy <- xy.coords(x, y, recycle=FALSE)
    N <- length(xy$x)
    # check for each polygon of polyregion if points are in the polygon
    locations <- sapply(polyregion@pts, function (poly) {
        pip <- point.in.polygon(xy$x, xy$y, poly$x, poly$y, mode.checked = mode.checked)
        if (poly$hole) { # if point is inside a hole then attribute -Inf
            ifelse(pip == 1, -Inf, 0)
        } else pip
    })
    if (N == 1) sum(locations) > 0 else
    .rowSums(locations, N, length(polyregion@pts)) > 0
}


### Maximum extent of a gpc.poly (i.e. maximum distance of two vertices)

maxExtent.gpc.poly <- function (object, ...)
{
    pts <- object@pts
    x <- unlist(lapply(pts, "[[", "x"), use.names=FALSE)
    y <- unlist(lapply(pts, "[[", "y"), use.names=FALSE)
    
    ## The diagonal of the bounding box provides a fast upper bound
    ##ext <- sqrt(diff(range(x))^2 + diff(range(y))^2)

    xy <- cbind(x,y)
    dists <- dist(xy)
    max(dists)
}


### Count number of instances at the same location of a SpatialPoint object

multiplicity.default <- function (x, ...)
{
    distmat <- as.matrix(dist(x))
    as.integer(rowSums(distmat == 0))
}

multiplicity.Spatial <- function (x, ...)
{
    multiplicity(coordinates(x))
}


### determine matrix with higher neighbourhood order based on spdep::nblag()
### given the binary matrix of first-order neighbours

nbOrder <- function (neighbourhood, maxlag = 1)
{
    if (!requireNamespace("spdep"))
        stop("package ", dQuote("spdep"),
             " is required to determine neighbourhood orders")

    stopifnot(isScalar(maxlag), maxlag > 0)
    checkNeighbourhood(neighbourhood)
    neighbourhood <- neighbourhood == 1           # convert to binary matrix
    nregions <- nrow(neighbourhood)
    maxlag <- as.integer(min(maxlag, nregions-1)) # upper bound of nb order
    
    if (maxlag == 1L) {
        storage.mode(neighbourhood) <- "integer"
        return(neighbourhood)
    }

    ## manually convert to spdep's "nb" class
    ## region.idxs <- seq_len(nregions)
    ## nb <- lapply(region.idxs, function(i) {
    ##     nbs <- which(neighbourhood[i,])
    ##     if (length(nbs) > 0L) nbs else 0L
    ## })
    ## class(nb) <- "nb"

    ## convert first-order neighbourhood to spdep's "nb" class
    nb <- spdep::mat2listw(neighbourhood)$neighbours
    attr(nb, "region.id") <- NULL

    ## compute higher order neighbours using spdep::nblag()
    nb.lags <- spdep::nblag(nb, maxlag=maxlag)

    ## Side note: fast method to determine neighbours _up to_ specific order:
    ## crossprod(neighbourhood) > 0  # up to second order neighbours (+set diag to 0)
    ## (neighbourhood %*% neighbourhood %*% neighbourhood) > 0  # up to order 3
    ## and so on...

    ## convert to a single matrix
    nbmat <- neighbourhood   # logical first-order matrix
    storage.mode(nbmat) <- "numeric"
    for (lag in 2:maxlag) {
        if (any(spdep::card(nb.lags[[lag]]) > 0L)) { # any neighbours of this order
            nbmat.lag <- spdep::nb2mat(nb.lags[[lag]], style="B",
                                       zero.policy=TRUE)
            nbmat <- nbmat + lag * nbmat.lag
        }
    }
    attr(nbmat, "call") <- NULL
    storage.mode(nbmat) <- "integer"

    ## message about maximum neighbour order by region
    maxlagbyrow <- apply(nbmat, 1, max)
    message("Note: range of maximum neighbour order by region is ",
            paste(range(maxlagbyrow), collapse="-"))

    ## Done
    nbmat
}

    
### determines which polygons of a SpatialPolygons object are at the border,
### i.e. have coordinates in common with the spatial union of all polygons

polyAtBorder <- function (SpP, snap=sqrt(.Machine$double.eps))
{
    W <- unionSpatialPolygons(SpP)      # length(W@polygons) == 1
    Wcoords <- unique(do.call("rbind",
                              lapply(W@polygons[[1]]@Polygons, coordinates)))
    atBorder <- sapply(SpP@polygons, function (x) {
        coords <- unique(do.call("rbind", lapply(x@Polygons, coordinates)))
        res <- FALSE
        for (i in seq_len(nrow(coords))) {
            if (any(spDistsN1(Wcoords, coords[i,]) < snap)) {
                res <- TRUE
                break
            }
        }
        res
    })
    names(atBorder) <- row.names(SpP)
    atBorder
}


### derive adjacency structure from SpatialPolygons
### (wrapping around functionality from the "spdep"-package)

poly2adjmat <- function (SpP, ..., zero.policy = TRUE)
{
    if (!requireNamespace("spdep"))
        stop("package ", dQuote("spdep"),
             " is required to derive adjacencies from SpatialPolygons")
    nb <- spdep::poly2nb(SpP, ...)
    adjmat <- spdep::nb2mat(nb, style="B", zero.policy=zero.policy)
    attr(adjmat, "call") <- NULL
    colnames(adjmat) <- rownames(adjmat)
    adjmat
}
