################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Auxiliary functions for operations on spatial data
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 666 $
### $Date: 2013-11-08 15:45:36 +0100 (Fre, 08 Nov 2013) $
################################################################################


### Polygonal Approximation of a Disc/Circle

discpoly <- function (center, radius, npoly = 64,
                      class = c("Polygon", "owin", "gpc.poly"),
                      hole = FALSE)
{
    class <- match.arg(class)
    if (class == "owin") { # use spatstat::disc
        res <- disc(radius=radius, centre=center, mask=FALSE, npoly=npoly)
        if (hole) {
            res$bdry[[1]]$x <- rev(res$bdry[[1]]$x)
            res$bdry[[1]]$y <- rev(res$bdry[[1]]$y)
            res$bdry[[1]]$hole <- TRUE
        }
        return(res)
    }

    ## do it myself for the "Polygon" and "gpc.poly" classes
    stopifnot(radius > 0, isScalar(npoly), npoly > 2)
    theta <- seq(2*pi, 0, length = npoly+1)[-(npoly+1)]   # for clockwise order
    if (hole) theta <- rev(theta)   # for anticlockwise order
    x <- center[1] + radius * cos(theta)
    y <- center[2] + radius * sin(theta)
    switch(class,
        "Polygon" = Polygon(cbind(c(x,x[1]),c(y,y[1])), hole=hole),
        "gpc.poly" = {
            pts <- list(list(x=x, y=y, hole=hole))
            if (isClass("gpc.poly") || requireNamespace("rgeos")) {
                new("gpc.poly", pts = pts)
            } else {
                warning("formal class \"gpc.poly\" not available")
                pts
            }
        }
    )
}


### Wrapper for polyclip or rgeos::gUnaryUnion or maptools::unionSpatialPolygons

unionSpatialPolygons <- function (SpP,
                                  method = c("rgeos", "polyclip", "gpclib"),
                                  ...)
{
    method <- match.arg(method)
    W <- switch(
        method,
        "polyclip" = {
            tiles_xylist <- xylist(SpP, reverse=FALSE)
            W_xylist <- polyclip::polyclip(tiles_xylist, tiles_xylist, "union",
                                           fillA = "nonzero", fillB = "nonzero",
                                           ...)
            ## FIXME: polyclip() seems to return owin-type vertex order?
            W_Polygons <- Polygons(
                lapply(W_xylist, function(p)
                       Polygon(cbind(p$x,p$y)[c(1L,length(p$x):1L),])),
                ID="1")
            SpatialPolygons(list(W_Polygons))
        },
        "rgeos" = rgeos::gUnaryUnion(SpP, ...),
        "gpclib" = {
            library("maptools") # loading only the namespace is not sufficient,
                                # since rgeosStatus is only set in .onAttach
            gpclibCheck() && maptools::gpclibPermit()
            maptools::unionSpatialPolygons(
                SpP, IDs = rep.int(1,length(SpP@polygons)),
                avoidGEOS = TRUE, ...)
        })
    ## ensure that W has exactly the same proj4string as SpP
    W@proj4string <- SpP@proj4string
    W
}

    
### Compute distance from points to boundary
### copied in part from function bdist.points() of the "spatstat" package
### authored by A. Baddeley and R. Turner (DEPENDS ON spatstat::distppl)
## xy is the coordinate matrix of the points
## poly is a polygonal domain of class "owin" or "gpc.poly"
## Note that we do not check if points are actually inside the polygonal domain
bdist <- function (xy, poly)
{
    result <- rep.int(Inf, nrow(xy))
    bdry <- if (is.polygonal(poly)) {
        poly$bdry
    } else if (inherits(poly, "gpc.poly")) {
        poly@pts
    } else stop("'poly' must be \"owin\" or \"gpc.poly\"")
    for (i in seq_along(bdry)) {
        polly <- bdry[[i]]
        px <- polly$x
        py <- polly$y
        nsegs <- length(px)
        for (j in seq_len(nsegs)) {
            j1 <- if (j < nsegs) j + 1L else 1L
            seg <- c(px[j], py[j], px[j1], py[j1])
            ## calculate distances of points to segment j of polygon i
            result <- pmin.int(result, distppl(xy, seg))
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

    
### determines which polygons of a SpatialPolygons object are at the border,
### i.e. have coordinates in common with the spatial union of all polygons

polyAtBorder <- function (SpP,
                          snap = sqrt(.Machine$double.eps),
                          method = "rgeos", ...)
{
    SpP <- as(SpP, "SpatialPolygons")
    W <- unionSpatialPolygons(SpP, method = method, ...)
    if (length(W@polygons) > 1)
        warning("unionSpatialPolygons() produced >1 Polygons-components")
    Wcoords <- unique(do.call("rbind",
                              lapply(W@polygons[[1]]@Polygons, coordinates)))
    atBorder <- sapply(SpP@polygons, function (x) {
        coords <- unique(do.call("rbind", lapply(x@Polygons, coordinates)))
        res <- FALSE
        for (i in seq_len(nrow(coords))) {
            if (any(spDistsN1(Wcoords, coords[i,], longlat=FALSE) < snap)) {
                res <- TRUE
                break
            }
        }
        res
    })
    names(atBorder) <- row.names(SpP)
    atBorder
}
