################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Wrappers of gpclib or rgeos functionality (polygons intersection/union),
### depending on surveillance.options("gpclib")
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 603 $
### $Date: 2013-07-18 11:54:49 +0200 (Don, 18 Jul 2013) $
################################################################################


### Compute the intersection of a "gpc.poly" or "SpatialPolygons" with a
### disc using gpclib or rgeos. The result is always a "gpc.poly".

intersectCircle <- function (Wgpc, W, center, r, npoly)
{
    gpclib <- gpclibCheck(fatal=FALSE)
    circle <- discpoly(center, r, npoly = npoly,
                       class = if (gpclib) "gpc.poly" else "Polygon")
    if (gpclib) {
        gpclib::intersect(circle, Wgpc)  # this order seems to be faster
    } else { # use rgeos
        circleSpP <- SpatialPolygons(list(Polygons(list(circle), "0")))
        circleSpP@proj4string <- W@proj4string
        res <- rgeos::gIntersection(circleSpP, W)
        as(res, "gpc.poly")      # coerce-method imported from rgeos via polyCub
    }
}


### Internal wrapper for maptools::unionSpatialPolygons

unionSpatialPolygons <- function (SpP)
{
    library("maptools") # loading only the namespace is not sufficient,
                        # since rgeosStatus is only set in .onAttach
    gpclib <- gpclibCheck(fatal=FALSE) && maptools::gpclibPermit()
    W <- maptools::unionSpatialPolygons(
        SpP, IDs = rep.int(1,length(SpP@polygons)),
        avoidGEOS = gpclib)
    ## ensure that W has exactly the same proj4string as SpP
    ## since the internal CRS()-call might have modified it
    W@proj4string <- SpP@proj4string
    W
}
