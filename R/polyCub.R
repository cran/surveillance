################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Methods for two-dimensional numerical integration over a polygonal domain
### cf. Section 3.2 of my Master's Thesis: http://epub.ub.uni-muenchen.de/11703/
###
### Copyright (C) 2009-2012 Sebastian Meyer
### $Revision: 463 $
### $Date: 2012-12-06 17:26:39 +0100 (Do, 06. Dez 2012) $
###
### The product Gauss cubature in 'polygauss()' is based on the corresponding
### MATLAB code ('polygauss') by Sommariva & Vianello (2007):
### "Product Gauss cubature over polygons based on Green's integration formula"
### Bit Numerical Mathematics, 47 (2), 441-453.
###
### PARAMS:
### polyregion: a polygon of class "gpc.poly" (gpclib::) or "owin" (spatstat::)
### f:  two-dimensional integrand. The function must take a coordinate matrix as
###     its first argument.
### ...: further arguments passed to f.
################################################################################


### Wrapper function for the various cubature methods
## ...: arguments of f or of the specific methods

polyCub <- function (polyregion, f, method = c("SV", "midpoint", "exact.Gauss"), ..., plot = FALSE)
{
	method <- match.arg(method)
	cl <- match.call()
	cl$method <- NULL
	cl[[1]] <- as.name(paste("polyCub", method, sep="."))
	if (method == "exact.Gauss") cl$f <- NULL
	int <- eval(cl, parent.frame())
	int  #structure(int, method = method)
}



########################################################################
### Alternative 1: Two-dimensional midpoint rule (using spatstat image)
########################################################################
## polyregion can be anything coercible to "owin" by as.owin
## eps: width and height of the pixels (squares)
## dimyx: number of subdivisions in each dimension

polyCub.midpoint <- function (polyregion, f, ..., eps = NULL, dimyx = NULL, plot = FALSE)
{
    # as.im needs seperate x and y arguments
    fxy <- function (x, y, ...) f(cbind(x,y), ...)

    # calculate pixel values of fxy
    IM <- tryCatch(
          as.im.function(X=fxy, W=polyregion, ..., eps=eps, dimyx=dimyx),
          error = function (e) {
              ## if eps was to small such that the dimensions of the image would
              ## be too big then the operation matrix(TRUE, nr, nc) throws an
              ## error. (try e.g. devnull <- matrix(TRUE, 1e6,1e6))
              ## unfortunately, it is not clear what we should do in this
              ## case... => stop
              stop("inapplicable choice of bandwidth (eps=", format(eps),
                   ") in midpoint rule:\n", e)
          })
    
    
### ILLUSTRATION ###
if (plot) {
    plot.im(IM, axes = TRUE, col=grey(31:4/35), main="")
    # add evaluation points (unsure about spatstat implementation of class "im")
    # both of the following commands worked with different versions of spatstat
    #with(IM, points(expand.grid(xcol, yrow), col=!is.na(v), cex=0.5))
    #with(IM, points(expand.grid(y=yrow, x=xcol)[2:1], col=!is.na(v), cex=0.5))
    plot(polyregion, add = TRUE, poly.args = list(lwd = 2))
}
####################
    
    # return the approximated integral
    pixelarea <- IM$xstep * IM$ystep
    int <- pixelarea * sum(IM$v, na.rm = TRUE)
    int
}



########################################################################
### Alternative 2: Product Gauss cubature of two-dimensional functions
### over simple polygons proposed by Sommariva & Vianello (2007)
########################################################################
## polyregion may be of classes "owin", "SpatialPolygons", or "gpc.poly"
## See function polygauss() below for further argument details

polyCub.SV <- function (polyregion, f, ...,
                        nGQ = 20, alpha = NULL, rotation = FALSE,
                        plot = FALSE)
{
    polys <- xylist(polyregion) # transform to something like "owin$bdry"
                                # which means anticlockwise vertex order with
                                # first vertex not repeated
    f <- match.fun(f)
    stopifnot(isScalar(nGQ), is.null(alpha) || (isScalar(alpha) && !is.na(alpha)))

    int1 <- function (poly) {
        nw <- polygauss(poly, nGQ, alpha, rotation)
        fvals <- f(nw$nodes, ...)
        cubature_val <- sum(nw$weights * fvals)
        ## if (!isTRUE(all.equal(0, cubature_val))) {
        ## if ((1 - 2 * as.numeric(poly$hole)) * sign(cubature_val) == -1)
        ## warning("wrong sign if positive integral")
        ## }
        cubature_val
    }
    respolys <- sapply(polys, int1, simplify = TRUE, USE.NAMES = FALSE)
    int <- sum(respolys)

    if (plot) {
        if (inherits(polyregion, "gpc.poly")) {
            gpclibCheck()
            plot(polyregion, poly.args=list(lwd=2), ann=FALSE)
        } else plot(polyregion, lwd=2, axes=TRUE, main="")
        for (i in seq_along(polys)) {
            nw <- polygauss(polys[[i]], nGQ, alpha, rotation)
            points(nw$nodes, cex=0.6, pch = i) #, col=1+(nw$weights<=0)
        }
    }

    int
}


### Function to calculate nodes and weights of the product Gauss cubature
### Code is based on the MATLAB implementation of Sommariva & Vianello (2007)
## Parameters:
## xy: list with elements "x" and "y" containing the polygon vertices in
##     _anticlockwise_ order (otherwise the result of the cubature will have a
##     negative sign) with first vertex not repeated at the end (like owin$bdry)
## N: degree of the one-dimensional Gauss quadrature rule
##    (see statmod::gauss.quad, on which this function depends)
## alpha: base-line at x = alpha (see the referenced paper for an explication).
##        If NULL (the default), the midpoint of the x-range is chosen if
##        rotation=FALSE and otherwise the x-coordinate of the rotated point P.
##        alpha = 0 seems to be reasonable if f has its maximum value at (0,0),
##        e.g. for the bivariate normal density with zero mean
## rotation: do efficiency rotation for convex polygons? (possibly a list of
##           points P and Q describing the preferred direction)
## kind: 1D quadrature rule (see statmod::gauss.quad)

polygauss <- function (xy, N, alpha = NULL, rotation = FALSE, kind = "legendre")
{
    ## convert to coordinate matrix
    xy <- cbind(x=xy[["x"]], y=xy[["y"]], deparse.level=0)

    
    ## POLYGON ROTATION
    
    xyrot <- if (identical(FALSE, rotation)) {
        if (is.null(alpha)) { # choose midpoint of x-range
            xrange <- range(xy[,1L])
            alpha <- (xrange[1L] + xrange[2L]) / 2
        }
        angle <- 0
        xy
    } else {
        if (identical(TRUE, rotation)) { # automatic choice of rotation angle
            ## such that for a convex polygon all nodes fall inside the polygon
            QP <- vertexpairmaxdist(xy)
            Q <- QP[1L,,drop=TRUE]
            P <- QP[2L,,drop=TRUE]
        } else if (is.list(rotation)) {  # predefined rotation
            stopifnot(is.vector(P <- rotation$P, mode="numeric") && length(P) == 2L,
                      is.vector(Q <- rotation$Q, mode="numeric") && length(Q) == 2L)
            stopifnot(any(P != Q))
            rotation <- TRUE
        } else {
            stop("'rotation' must be logical or a list of points \"P\" and \"Q\"")
        }
        rotmat <- rotmatPQ(P,Q)
        angle <- attr(rotmat, "angle")
        if (is.null(alpha)) {
            Prot <- rotmat %*% P
            alpha <- Prot[1]
        }
        xy %*% t(rotmat)   # = t(rotmat %*% t(xy))
    }

    
    ## COMPUTE NODES AND WEIGHTS OF 1D GAUSS QUADRATURE RULE.
    
    ## DEGREE "N" (as requested) (ORDER GAUSS PRIMITIVE)
    nw_N <- statmod::gauss.quad(n = N, kind = kind)

    ## DEGREE "M" = N+1 (ORDER GAUSS INTEGRATION)
    M <- N + 1L
    nw_M <- statmod::gauss.quad(n = M, kind = kind)

    
    ## COMPUTE 2D NODES AND WEIGHTS.

    xbd <- xyrot[,1L,drop=TRUE]
    ybd <- xyrot[,2L,drop=TRUE]
    nw <- .polygauss(c(xbd,xbd[1L]), c(ybd,ybd[1L]), alpha,
                     nw_N$nodes, nw_N$weights, nw_M$nodes, nw_M$weights)
    #nw <- .Call("polygauss",
    #            xyrot[,1L],xyrot[,2L],alpha,
    #            nw_N$nodes,nw_N$weights,nw_M$nodes,nw_M$weights,
    #            PACKAGE="surveillance")

    ## back-transform rotated nodes by t(t(rotmat) %*% t(nodes))
    ## (inverse of rotation matrix is its transpose)
    nodes <- if (rotation) nw$nodes %*% rotmat else nw$nodes
    
    ## Done.
    list(nodes=nodes, weights=nw$weights, angle = angle, alpha = alpha)
}


### Main part which calculates nodes and weights of the product Gauss cubature
### Code is based on the MATLAB implementation of Sommariva & Vianello (2007)
## TODO: The efficient implementation of this function in C
## would increase the speed of the cubature
## Parameters:
## x, y: coordinates of the polygon's vertices in _anticlockwise_ order
##       (otherwise the result of the cubature will have a negative sign)
##       with _REPEATED FIRST VERTEX_ at the end
## alpha: "base-line"
## s/w_N/M: nodes and weights of univariate Gauss rules of orders N and M

.polygauss <- function (x, y, alpha, s_N, w_N, s_M, w_M)
{
    L <- length(x) - 1L                 # number of sides of the polygon
    N <- length(s_N)
    M <- length(s_M)
    maxNodes <- L*M*N
    nodes <- matrix(0, nrow=maxNodes, ncol=2L)
    weights <- numeric(maxNodes)
    K <- 0L

    for (side in seq_len(L))
    {
        x1 <- x[side];  x2 <- x[side+1L]
        y1 <- y[side];  y2 <- y[side+1L]

        if ((x1 == alpha && x2 == alpha) || (y2 == y1)) {
            ## skip: side lies on base-line or is orthogonal to it
            next
        }

        if (x2 == x1) { # side is parallel to base-line => degree N
            degree_loc <- N
            s_loc <- s_N
            w_loc <- w_N
        } else { # degree M=N+1
            degree_loc <- M
            s_loc <- s_M
            w_loc <- w_M
        }

        nw_loc <- .polygauss.side(x1, y1, x2, y2, s_loc, w_loc, s_N, w_N, alpha)

        ## add nodes and weights for this side of the polygon
        nNodes_loc <- degree_loc * N
        rowidx <- K + seq_len(nNodes_loc)
        nodes[rowidx,] <- nw_loc$nodes
        weights[rowidx] <- nw_loc$weights
        K <- K + nNodes_loc
    }
        
    # only the first K entries of 'weights' and rows of 'nodes_x' and 'nodes_y'
    # have been filled, the remainder until 'maxRows', contains the zeros
    # from initialisation.
    seqK <- seq_len(K)
    weights <- weights[seqK]
    nodes <- nodes[seqK,,drop=FALSE]

    ## Done
    ret <- list(nodes = nodes, weights = weights)
    ret
}

.polygauss.side <- function (x1, y1, x2, y2, s_loc, w_loc, s_N, w_N, alpha)
{
    degree_loc <- length(s_loc)
    N <- length(s_N)
    
    half_pt_x <- (x1+x2)/2;  half_length_x <- (x2-x1)/2;
    half_pt_y <- (y1+y2)/2;  half_length_y <- (y2-y1)/2;
    
    ## GAUSSIAN POINTS ON THE SIDE.
    x_gauss_side <- half_pt_x + half_length_x * s_loc
    y_gauss_side <- half_pt_y + half_length_y * s_loc

    ## construct weights
    scaling_fact_minus <- (x_gauss_side - alpha) / 2
    weights1 <- (half_length_y * scaling_fact_minus) * w_loc # length=degree_loc
    weights <- tcrossprod(weights1, w_N) # degree_loc x N
    #weights <- rep(w_N, each=degree_loc) * rep.int(weights1, N)  # lenth=degree_loc*N
    
    ## construct nodes: x and y coordinates ARE STORED IN MATRICES.
    ## A COUPLE WITH THE SAME INDEX IS A POINT, i.e. P_i=(x(k),y(k)).
    term_2 <- matrix(scaling_fact_minus, nrow = degree_loc, ncol = N)
    rep_n_Np1 <- matrix(s_N + 1, nrow = degree_loc, ncol = N, byrow = TRUE)
    nodes_x <- rep_n_Np1 * term_2 + alpha # degree_loc x N
    #nodes_y <- matrix(y_gauss_side, nrow = degree_loc, ncol = N)
    nodes_y <- rep.int(y_gauss_side, N)   # no matrix dimensions
    
    ## Done (return as 2-column matrix, and weights vector)
    list(nodes = cbind(x=c(nodes_x), y=nodes_y), weights = c(weights))
}

vertexpairmaxdist <- function (xy)
{
    ## compute euclidean distance matrix
    distances <- dist(xy)
    size <- attr(distances, "Size")
    
    ## select two points with maximum distance
    maxdistidx <- which.max(distances)
    lowertri <- seq_along(distances) == maxdistidx
    mat <- matrix(FALSE, size, size)
    mat[lower.tri(mat)] <- lowertri
    QPidx <- which(mat, arr.ind=TRUE, useNames=FALSE)[1L,]
    xy[QPidx,]    
}

rotmatPQ <- function (P, Q)
{
    direction_axis <- (Q-P) / sqrt(sum((Q-P)^2))
    
    ## determine rotation angle
    rot_angle_x <- acos(direction_axis[1L])
    rot_angle_y <- acos(direction_axis[2L])
    
    rot_angle <- if (rot_angle_y <= pi/2) {
        if (rot_angle_x <= pi/2) -rot_angle_y else rot_angle_y
    } else {
        if (rot_angle_x <= pi/2) pi-rot_angle_y else rot_angle_y
    }
    ## cat(sprintf(' [ANGLE CLOCKWISE (IN DEGREES)]: %5.5f\n', rot_angle*180/pi))

    ## rotation matrix
    rot_matrix <- diag(cos(rot_angle), nrow=2L)
    rot_matrix[2:3] <- c(-1,1) * sin(rot_angle) # clockwise rotation
    structure(rot_matrix, angle=rot_angle)
}



#######################################################################
### Alternative 3: Quasi-exact cubature of the bivariate normal density
### based on triangulation and formulae from Chapter 26 of the famous
### Abramowitz & Stegun handbook (Section 26.9, Example 9, pp. 956f.)
#######################################################################
## this is quite complicated because A&S formula is only for triangles where one
## vertex is the origin (0,0)
## For each triangle of the tristrip we have to check in which of the 6 outer regions
## of the triangle the origin (0,0) lies and adapt the signs in the formula adequately
## (AOB+BOC-AOC) or (AOB-AOC-BOC) or (AOB+AOC-BOC) or (AOC+BOC-AOB) or ...
## In contrast, the most time consuming step is the evaluation of pmvnorm in .V

## helper function, which calculates the integral of the standard bivariat normal
## over a triangle bounded by y=0, y=ax, x=h (cf. formula 26.3.23)
.V <- function(h,k) {
    a <- k/h
    rho <- -a/sqrt(1+a^2)
    # V = 0.25 + L(h,0,rho) - L(0,0,rho) - Q(h) / 2
    # L(0,0,rho) = 0.25 + asin(rho) / (2*pi)
    # V = L(h,0,rho) - asin(rho)/(2*pi) - Q(h) / 2
    Lh0rho <- mvtnorm::pmvnorm(
        lower = c(h,0), upper = c(Inf,Inf), mean = c(0,0), corr = matrix(c(1,rho,rho,1),2,2)
    )
    Qh <- pnorm(h, mean = 0, sd = 1, lower.tail = FALSE)
    return(Lh0rho - asin(rho)/2/pi - Qh/2)
}

## helper function, which calculates the integral of the standard bivariat normal over a triangle A0B
.intTriangleAS0 <- function (A, B)
{
    d <- sqrt(sum((B-A)^2))
    h <- abs(B[2]*A[1] - A[2]*B[1]) / d
    if (h == 0) return(0)
    k1 <- abs(A[1]*(B[1]-A[1]) + A[2]*(B[2]-A[2])) / d
    k2 <- abs(B[1]*(B[1]-A[1]) + B[2]*(B[2]-A[2])) / d
    
    V2 <- .V(h, k2)
    V1 <- .V(h, k1)
    res <- if (isTRUE(all.equal(k1+k2, d))) V2 + V1
        else if (isTRUE(all.equal(abs(k2-k1), d))) abs(V2 - V1)
        else stop("something went wrong...")
    attr(res, "error") <- attr(V1, "error") + attr(V2, "error")
    return(res)
}

## checks if point1 and point2 lie on the same side of a line through linepoint1 and linepoint2
## see, e.g., http://www.gamedev.net/community/forums/topic.asp?topic_id=457450
.pointsOnSameSide <- function (linepoint1, linepoint2, point1, point2 = c(0,0))
{
    n <- c(-1,1) * rev(linepoint2-linepoint1)   # normal vector
    S <- dotprod(point1-linepoint1,n) * dotprod(point2-linepoint1,n)
    return(S > 0)
}

## helper function, which calculates the integral of the standard bivariat normal over a triangle ABC
.intTriangleAS <- function (xy)
{
    A <- xy[1,]
    B <- xy[2,]
    C <- xy[3,]
    intAOB <- .intTriangleAS0(A, B)
    intBOC <- .intTriangleAS0(B, C)
    intAOC <- .intTriangleAS0(A, C)
    
    # determine signs of integrals
    signAOB <- -1 + 2*.pointsOnSameSide(A,B,C)
    signBOC <- -1 + 2*.pointsOnSameSide(B,C,A)
    signAOC <- -1 + 2*.pointsOnSameSide(A,C,B)
    
    int <- signAOB*intAOB + signBOC*intBOC + signAOC*intAOC
    attr(int, "error") <- attr(intAOB, "error") + attr(intBOC, "error") + attr(intAOC, "error")
    return(int)
}

## main function to be called by the user
polyCub.exact.Gauss <- function (polyregion, mean = c(0,0), Sigma = diag(2), plot = FALSE)
{
    gpclibCheck()
    polyregion <- as(polyregion, "gpc.poly")
    
    # coordinate transformation so that the standard bivariat normal density
    # can be used in integrations (cf. formula 26.3.22)
    rho <- cov2cor(Sigma)[1,2]
    sdx <- sqrt(Sigma[1,1])
    sdy <- sqrt(Sigma[2,2])
    polyregion@pts <- lapply(polyregion@pts, function (poly) {
        list(x = ((poly$x-mean[1])/sdx + (poly$y-mean[2])/sdy) / sqrt(2+2*rho),
             y = ((poly$y-mean[2])/sdy - (poly$x-mean[1])/sdx) / sqrt(2-2*rho),
             hole = poly$hole)
    })
    
    # triangulation: tristrip returns a list where each element is a coordinate matrix of vertices of triangles
    triangleSets <- gpclib::tristrip(polyregion)
    
### ILLUSTRATION ###
if (plot) {
    plot(polyregion, poly.args=list(lwd=2), ann=FALSE)
    lapply(triangleSets, lines, lty=2)
}
####################

    integrals <- sapply(triangleSets, function (triangles) {
        int <- 0
        error <- 0
        nTriangles <- nrow(triangles) - 2
        for (i in seq_len(nTriangles)) {
            res <- .intTriangleAS(triangles[i+(0:2),])
            err <- attr(res, "error")
            int <- int + res
            if (length(err) == 1L) error <- error + err   # sometimes err==numeric(0) (probably meaning err=0)
        }
        c(int, nTriangles, error)
    })
    int <- sum(integrals[1,])
    # number of .V() evaluations (if 'h' in .intTriangleAS0 was always different from 0)
    attr(int, "nEval") <- 6 * sum(integrals[2,])
    # approximate absolute integration error
    attr(int, "error") <- sum(integrals[3,])
    return(int)
}

