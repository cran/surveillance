################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Auxiliary functions for operations on spatial data
###
### Copyright (C) 2009-2014 Sebastian Meyer
### $Revision: 1086 $
### $Date: 2014-10-23 20:33:18 +0200 (Thu, 23 Oct 2014) $
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


### sample n points uniformly on a disc with radius r

runifdisc <- function (n, r = 1, buffer = 0)
{
    stopifnot(buffer <= r)
    rangle <- runif(n, 0, 2*pi)
    rdist <- r * sqrt(runif(n, (buffer/r)^2, 1))
    rdist * cbind(cos(rangle), sin(rangle))
}


### Count number of instances at the same location of a SpatialPoints object
## NOTE: the default multiplicity-method has been integrated into the spatstat
## package which we import

multiplicity.Spatial <- function (x) multiplicity(coordinates(x))

    
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


### sp.layout item for spplot() to draw labels for Spatial* objects

layout.labels <- function (obj, labels = TRUE)
{
    stopifnot(inherits(obj, "Spatial"))

    ## get region labels
    getLabels <- function (labels) {
        if (isTRUE(labels)) {
            row.names(obj)
        } else if (length(labels) == 1L &&
                   (is.numeric(labels) | is.character(labels))) {
            if (!"data" %in% slotNames(obj))
                stop("no data slot to select labels from")
            obj@data[[labels]]
        } else labels
    }
    
    ## convert labels argument to a list
    labels.args <- if (is.list(labels)) {
        labels
    } else if (!is.null(labels) && !identical(labels, FALSE)) {
        list(labels = getLabels(labels))
    } else { # labels = FALSE or labels = NULL
        return(NULL)
    }

    ## set default coordinates for panel.text() and parse labels
    labels.args <- modifyList(list(x = coordinates(obj), labels = TRUE),
                              labels.args)
    labels.args$labels <- getLabels(labels.args$labels)

    ## return layout item
    c("panel.text", labels.args)
}
