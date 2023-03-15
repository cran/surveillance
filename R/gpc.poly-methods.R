################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Methods for gpc.poly polygons
### These are no longer used by the surveillance package itself
###
### Copyright (C) 2009-2013 Sebastian Meyer
### $Revision: 2946 $
### $Date: 2023-03-08 21:52:14 +0100 (Wed, 08. Mar 2023) $
################################################################################


### Redefinition of gpclib's scale.poly method to also do centering

scale.gpc.poly <-
    function (x, center = c(0,0), scale = c(1,1)) {
        gpcWarning()
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
    gpcWarning()
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
