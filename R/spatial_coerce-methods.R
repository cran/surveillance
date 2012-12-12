################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Additional coerce-methods between different spatial classes
###
### Copyright (C) 2012 Sebastian Meyer
### $Revision: 467 $
### $Date: 2012-12-10 16:36:54 +0100 (Mo, 10. Dez 2012) $
################################################################################

### Note the varying polygon specifications in the packages:
## sp: REPEAT first vertex at the end (closed)
##     anticlockwise = hole, clockwise = normal boundary
## spatstat: do NOT REPEAT first vertex
##     anticlockwise = normal boundary, clockwise = hole
## gpc.poly: NOT DOCUMENTED, but it seems that it is prefered to not repeat the
##     first vertex and to have clockwise vertex order for the normal boundary.
##     The coerce-methods from "[Spatial]Polygons" or "owin" to "gpc.poly"
##     defined below follow this convention for "gpc.poly". 



### xylist methods get the simple list of polygons from the various classes
### where each component (polygon) is a simple list of "x" and "y" coordinates,
### which give the coordinates of the vertices of the polygon
### following the "owin" convention (anticlockwise order without repeating any
### vertex). There may be additional elements "area" and "hole" in each
### component, but these are not necessary.

xylist.owin <- function (object, ...) object$bdry
xylist.gpc.poly <- function (object, ...) xylist.owin(gpc2owin(object))
xylist.SpatialPolygons <- function (object, ...)
    xylist.owin(maptools::as.owin.SpatialPolygons(object))
## for the default method, no transformation is performed, we only check that
## polys are not closed (first vertex not repeated)
xylist.default <- function (object, ...) {
    lapply(object, function(xy) {
        poly <- xy.coords(xy)[c("x","y")]
        if (isClosed(poly)) {
            n <- length(poly$x)
            sel <- seq_len(n-1L)
            poly$x <- poly$x[sel]
            poly$y <- poly$y[sel]
        }
        poly
    })
}


### Method for coercion from "Polygons" (sp) to "gpc.poly" (gpclib) and vice versa

setAs(from = "Polygons", to = "gpc.poly", def = function (from)
    {
        gpclibCheck()
        pls <- slot(from, "Polygons")
        pts <- lapply(pls, function (sr) {
            coords <- coordinates(sr)
            n <- nrow(coords) - 1   # number of vertices
            list(x = coords[seq_len(n),1],
                 y = coords[seq_len(n),2],
                 hole = sr@hole)
        })
        new("gpc.poly", pts = pts)
    }
)

setAs(from = "gpc.poly", to = "Polygons", def = function (from)
    {
        srl <- lapply(from@pts, function (poly) {
            if (isClosed(poly)) {
                Polygon(cbind(poly$x,poly$y), hole = poly$hole)
            } else {
                Polygon(cbind(c(poly$x,poly$x[1]), c(poly$y,poly$y[1])),
                        hole = poly$hole)
            }
        })
        Polygons(srl, ID = "1")
    }
)



### Method for coercion from "SpatialPolygons" (sp) to "gpc.poly" (gpclib)
### This method also applies to "SpatialPolygonsDataFrame" (inherited)

setAs(from = "SpatialPolygons", to = "gpc.poly", def = function (from)
    {
        gpclibCheck()
        polygonsList <- polygons(from)@polygons
        gpc <- new("gpc.poly")
        for (i in seq_along(polygonsList))
        {
            gpc <- gpclib::append.poly(gpc, as(polygonsList[[i]], "gpc.poly"))
        }
        gpc
    }
)



### Method for coercion from "owin" (spatstat) to "gpc.poly" (gpclib)
### Coercion the other way round is defined in spatstat::gpc2owin

setAs(from = "owin", to = "gpc.poly", def = function (from)
    {
        gpclibCheck()
        pts <- lapply(from$bdry, function (poly) {
            list(x = rev(poly$x), y = rev(poly$y), hole = poly$hole)
        })
        new("gpc.poly", pts = pts)
    }
)

