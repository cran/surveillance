################################################################################
### Auxiliary functions for operations on spatial data
###
### Copyright (C) 2009-2015,2018,2021-2024  Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################


### Polygonal Approximation of a Disc/Circle

discpoly <- function (center, radius, npoly = 64,
                      class = c("Polygon", "owin", "gpc.poly"),
                      ## FIXME: should add polygonal "sfg" (or "sfc"?)
                      hole = FALSE)
{
    class <- match.arg(class)
    if (class == "owin") { # use spatstat.geom::disc
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
    theta <- seq(2*pi, 0, length.out = npoly+1)[-(npoly+1)]   # for clockwise order
    if (hole) theta <- rev(theta)   # for anticlockwise order
    x <- center[1] + radius * cos(theta)
    y <- center[2] + radius * sin(theta)
    switch(class,
        "Polygon" = Polygon(cbind(c(x,x[1]),c(y,y[1])), hole=hole),
        "gpc.poly" = {
            pts <- list(list(x=x, y=y, hole=hole))
            if (isClass("gpc.poly")) { #|| requireNamespace("gpclib")
                new("gpc.poly", pts = pts)
            } else {
                warning("formal class \"gpc.poly\" not available")
                pts
            }
        }
    )
}


### Wrapper for polyclip or sf::st_union

unionSpatialPolygons <- function (SpP,
                                  method = c("sf", "polyclip"),
                                  ...)
{
    if (identical(method, "gpclib") || identical(method, "rgeos")) {
        .Deprecated(msg = sprintf("method = \"%s\" is retired; using default", method))
        method <- NULL
    }
    method <- match.arg(method)
    W <- switch(
        method,
        "sf" = {
            sf::as_Spatial(sf::st_union(sf::st_as_sfc(SpP), ...),
                           IDs = "1") # as in rgeos::gUnaryUnion() before
        },
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
        }
    )
    ## ensure that W has exactly the same proj4string as SpP
    W@proj4string <- SpP@proj4string
    W
}


### internal implementation of as(W, "owin") from polyCub
### to avoid upgrade problems with polyCub <= 0.7.1 referring to old spatstat
### and to avoid calling as(W, "owin") from maptools (depends on load order)

SpP2owin <- function (W, ...) owin(poly = xylist(W), ...)


### Compute distance from points to a polygonal boundary

## since spatstat 1.56-0, bdist.points() interfaces C-code via
## spatstat.utils:::distppllmin, which is faster than nncross.ppp()
bdist <- function (xy, poly)  # poly is a polygonal "owin"
{
    bdist.points(ppp(x = xy[,1L], y = xy[,2L], window = poly, check = FALSE))
}
## Example: bdist(coordinates(imdepi$events), as(imdepi$W, "owin"))


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
                          method = "sf", ...)
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


### sp.layout items for spplot()

## draw labels for Spatial* objects
layout.labels <- function (obj, labels = TRUE, plot = FALSE)
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

    if (plot) {
        ## plot labels in the traditional graphics system
        do.call("text", labels.args)
    } else {
        ## return layout item for use by spplot()
        c("panel.text", labels.args)
    }
}

## draw a scalebar with labels
layout.scalebar <- function (obj, corner = c(0.05, 0.95), scale = 1,
                             labels = c(0, scale), height = 0.05,
                             pos = 3, ..., plot = FALSE)
{
    stopifnot(inherits(obj, "Spatial"))
    BB <- bbox(obj)
    force(labels)  # the default should use the original 'scale' value in km
    if (isFALSE(is.projected(obj, warn = TRUE))) {
        ## 'obj' has longlat coordinates, 'scale' is interpreted in kilometres
        scale <- .scale2longlat(t(rowMeans(BB)), scale)
    }
    offset <- BB[, 1L] + corner * apply(BB, 1L, diff.default)
    textfun <- if (plot) "text" else "panel.text"
    lis <- list(
        list("SpatialPolygonsRescale", layout.scale.bar(height = height),
             offset = offset, scale = scale, fill = c(NA, 1),
             plot.grid = !plot),
        list(textfun, x = offset[1L], y = offset[2L],
             labels = labels[1L], pos = pos, ...),
        list(textfun, x = offset[1L] + scale[1L], y = offset[2L],
             labels = labels[2L], pos = pos, ...)
    )
    if (plot) {
        for (li in lis) eval(do.call("call", li))
    } else {
        lis
    }
}

.scale2longlat <- function (focusLL, distKM)
{
    ## .destPoint() is copied from the "raster" package by Robert J. Hijmans
    ## 'p' is a longlat coordinate matrix, 'd' is a vector of distances in metres
    .destPoint <- function (p, d, b=90, r=6378137) {
        toRad <- pi/180
        lon1 <- p[, 1] * toRad
        lat1 <- p[, 2] * toRad
        b <- b * toRad
        lat2 <- asin(sin(lat1) * cos(d/r) + cos(lat1) * sin(d/r) * cos(b))
        lon2 <- lon1 + atan2(sin(b) * sin(d/r) * cos(lat1), cos(d/r) - sin(lat1) * sin(lat2))
        lon2 <- (lon2 + pi)%%(2 * pi) - pi
        cbind(lon2, lat2)/toRad
    }
    rightLL <- .destPoint(focusLL, distKM * 1000)
    rightLL[,1L] - focusLL[,1L]
}

## internal wrapper for sp::is.projected to catch its sf dependence
is.projected <- function (obj, warn = FALSE)
{
    res <- tryCatch(sp::is.projected(obj), packageNotFoundError = identity)
    ## should no longer be needed as the sp version has kept the no-sf fallback:
    if (inherits(res, "packageNotFoundError")) {
        ## fallback: grep for longlat in (deprecated) Proj.4 representation
        p4s <- as.character(obj@proj4string@projargs)
        res <- if (is.na(p4s) || !nzchar(p4s)) NA
               else !grepl("longlat", p4s, fixed = TRUE)
    }
    if (is.na(res) && warn)
        warning("could not determine projection status")
    res
}

## internal replacement for sp::mapasp using the above is.projected wrapper
mapasp <- function (data)
{
    if (isFALSE(is.projected(data))) {
        bb <- bbox(data)
        xlim <- bb[1L,]
        ylim <- bb[2L,]
        (diff(ylim)/diff(xlim)) / cos((mean(ylim) * pi)/180)
    } else "iso"
}


### determine the total area of a SpatialPolygons object

areaSpatialPolygons <- function (obj, byid = FALSE)
{
    ## if (requireNamespace("rgeos", quietly = TRUE)) {
    ##     rgeos::gArea(obj, byid = byid)
    ## } else {
        areas <- vapply(
            X = obj@polygons,
            FUN = function (p) sum(
                vapply(X = p@Polygons,
                       FUN = function (x) (1-2*x@hole) * x@area,
                       FUN.VALUE = 0, USE.NAMES = FALSE)
            ),
            FUN.VALUE = 0, USE.NAMES = FALSE
        )
        if (byid) setNames(areas, row.names(obj)) else sum(areas)
    ## }
}


### build a SpatialGrid of approx. 'n' square cells covering 'bbox(obj)'
## replacement for maptools::Sobj_SpatialGrid()

makeGrid <- function (obj, n = 128)
{
    bb <- bbox(obj)
    xlim <- bb[1L,]
    ylim <- bb[2L,]
    xr <- xlim[[2L]] - xlim[[1L]]
    yr <- ylim[[2L]] - ylim[[1L]]
    asp <- if (lscape <- xr > yr) xr/yr else yr/xr
    n1 <- ceiling(sqrt(n*asp))
    n2 <- ceiling(n1/asp)
    size <- (if (lscape) xr else yr) / n1
    offset <- bb[,1L] + size/2
    grd <- GridTopology(offset, c(size, size),
                        if (lscape) c(n1, n2) else c(n2, n1))
    SpatialGrid(grd, proj4string = obj@proj4string)
}

if (FALSE) {
## alternative that usually does not produce square cells,
## but has the same extent as bbox(obj),
## similar to sf::st_make_grid(obj, n = magic.dim(.));
## this would produce different spatial intensity plots than before
makeGrid <- function (obj, n = 128)
{
    bb <- bbox(obj)
    xlim <- bb[1L,]
    ylim <- bb[2L,]
    xr <- xlim[[2L]] - xlim[[1L]]
    yr <- ylim[[2L]] - ylim[[1L]]
    if (length(n) == 1) {
        n <- magic.dim(n)
        if (xr > yr) n <- rev(n)
    } else stopifnot(length(n) == 2)
    size <- c(xr, yr)/n
    offset <- bb[, 1L] + size/2
    grd <- GridTopology(offset, size, n)
    SpatialGrid(grd, proj4string = obj@proj4string)
}
}
