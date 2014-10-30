################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Compute distance from points to a polygonal boundary
###
### Copyright (C) 2009-2014 Sebastian Meyer
### The functions bdist() and distppl() are adapted from
### spatstat::bdist.points and spatstat:::distppl authored by Adrian Baddeley,
### who again took distppl() from an original by Rob Foxall, 1997
###
### $Revision: 1086 $
### $Date: 2014-10-23 20:33:18 +0200 (Thu, 23 Oct 2014) $
################################################################################


## xy is the coordinate _matrix_ of the points
## poly is a polygonal domain of class "owin" or "gpc.poly"
## Note that we do not check if points are actually inside the polygonal domain
bdist <- function (xy, poly)
{
    result <- rep.int(Inf, length(xy)/2) # faster than nrow, xy must be a matrix
    bdry <- if (is.polygonal(poly)) {
        poly$bdry
    } else if (inherits(poly, "gpc.poly")) {
        poly@pts
    } else stop("'poly' must be a polygonal \"owin\" or a \"gpc.poly\"")
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

## compute distances from each of a list of points p[i,]
## to a single line segment l = (x1, y1, x2, y2)
distppl <- function (p, l)
{
    xp <- p[,1L]
    yp <- p[,2L]
    dx <- l[3L] - l[1L]
    dy <- l[4L] - l[2L]
    leng <- sqrt(dx^2 + dy^2)
    ## vector from 1st endpoint to p
    xpl <- xp - l[1L]
    ypl <- yp - l[2L]
    ## distance from p to 1st & 2nd endpoints
    d1 <- sqrt(xpl^2 + ypl^2)
    d2 <- sqrt((xp - l[3L])^2 + (yp - l[4L])^2)
    dmin <- pmin.int(d1, d2)
    ## test for zero length
    if(leng < .Machine$double.eps)
        return(dmin)
    ## rotation sine & cosine
    co <- dx/leng
    si <- dy/leng
    ## back-rotated coords of p
    xpr <- co * xpl + si * ypl
    ypr <-  - si * xpl + co * ypl
    ## ypr is perpendicular distance to infinite line
    ## Applies only when xp, yp in the middle
    middle <- (xpr >= 0 & xpr <= leng)
    if (any(middle))
        dmin[middle] <- pmin.int(dmin[middle], abs(ypr[middle]))
    
    return(dmin)
}
