################################################################################
### Spatial helper functions
###
### Author: Sebastian Meyer
### $Date: 2012-07-23 17:58:36 +0200 (Mo, 23. Jul 2012) $
################################################################################


## Returns a Polygon representing a disc (in planar coordinates)

# center: center of the disc
# r: radius in km
# npoly: Number of edges of the polygonal approximation
# hole: hole flag of the polygon
discpoly <- function (center, r, npoly = 64,
    class = c("gpc.poly", "owin", "Polygon"), hole = FALSE)
{
    class <- match.arg(class)
    if (class == "owin") {
        res <- spatstat::disc(radius = r, centre = center, mask = FALSE, npoly = npoly)
        if (hole) {
            res$bdry[[1]]$x <- rev(res$bdry[[1]]$x)
            res$bdry[[1]]$y <- rev(res$bdry[[1]]$y)
            res$bdry[[1]]$hole <- TRUE
        }
        return(res)
    }

    theta <- seq(2*pi, 0, length = npoly+1)[-(npoly+1)]   # for clockwise order
    if (hole) theta <- rev(theta)   # for anticlockwise order
    x <- center[1] + r * cos(theta)
    y <- center[2] + r * sin(theta)
    switch(class,
        "Polygon" = Polygon(cbind(c(x,x[1]),c(y,y[1])), hole=hole),
        "gpc.poly" = new("gpc.poly", pts = list(list(x=x, y=y, hole=hole)))
    )
}


### Redefinition of the gpclib::scale.poly method for gpc.poly's to also do centering

setMethod("scale.poly", signature(x = "gpc.poly"),
# alternatively: define S3 method for scale
#scale.gpc.poly <-
    function (x, center = c(0,0), scale = c(1,1)) {
        x@pts <- lapply(x@pts, function (p) {
            p$x <- (p$x-center[1]) / scale[1]
            p$y <- (p$y-center[2]) / scale[2]
            p
        })
        x
    }
)


### Same as inside.owin for gpc.poly (using point.in.polygon from package sp)

inside.gpc.poly <- function(x, y = NULL, polyregion, mode.checked = FALSE)
{
	xy <- xy.coords(x, y, recycle=FALSE)
	N <- length(xy$x)
    # check for each polygon of polyregion if points are in the polygon
    locations <- sapply(gpclib::get.pts(polyregion), function (poly) {
        pip <- sp::point.in.polygon(xy$x, xy$y, poly$x, poly$y, mode.checked = mode.checked)
        if (poly$hole) { # if point is inside a hole then attribute -Inf
            ifelse(pip == 1, -Inf, 0)
        } else pip
    })
    inside <- if (N == 1L) sum(locations) > 0 else rowSums(locations) > 0
    return(inside)
}


## Count number of instances at the same location of a SpatialPoint object

multiplicity.default <- function (x, ...)
{
    distmat <- as.matrix(dist(x))
    as.integer(rowSums(distmat == 0))
}

multiplicity.Spatial <- function (x, ...)
{
    multiplicity(coordinates(x))
}


### Checks if the first and last coordinates of a coordinate matrix are equal

isClosed <- function (coords)
{
    xycoords <- xy.coords(coords)[c("x","y")]
    n <- length(xycoords$x)
    return(identical(xycoords$x[1], xycoords$x[n]) &&
           identical(xycoords$y[1], xycoords$y[n]))
}


## Compute the intersection of a "gpc.poly" with a "discpoly" (and centering)

intersectCircle <- function (Wgpc, center, r, npoly)
{
    circle <- discpoly(center = center, r = r, npoly = npoly,
                       class = "gpc.poly", hole = FALSE)
    intersection <- gpclib::intersect(circle, Wgpc)  # this order seems to be faster
    scale.poly(intersection, center = center) # use scale method as defined above
}

