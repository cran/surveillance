################################################################################
### coerce-methods for the different spatial classes
###
### Author: Sebastian Meyer
### $Date: 2012-06-27 15:25:47 +0200 (Mi, 27. Jun 2012) $
################################################################################

### Note the varying polygon specifications in the packages:
# sp: REPEAT first vertex at the end (closed)
#     anticlockwise = hole, clockwise = normal boundary
# spatstat: do NOT REPEAT first vertex
#     anticlockwise = normal boundary, clockwise = hole
# gpc.poly: NOT DOCUMENTED, but it seems that it is prefered to not repeat the
#     first vertex and to have clockwise vertex order for the normal boundary. 


### Method for coercion from "Polygons" (sp) to "gpc.poly" (gpclib) and vice versa

setAs(from = "Polygons", to = "gpc.poly", def = function (from)
    {
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
        srl <- lapply(gpclib::get.pts(from), function (poly) {
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
### Coercion the other way round is defined in spatstat::as.owin.gpc.poly

setAs(from = "owin", to = "gpc.poly", def = function (from)
    {
        pts <- lapply(from$bdry, function (poly) {
            list(x = rev(poly$x), y = rev(poly$y), hole = poly$hole)
        })
        new("gpc.poly", pts = pts)
    }
)

