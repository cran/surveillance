################################################################################
### Methods for two-dimensional numerical integration over a polygonal domain
### (see Section 3.2 of my Master's Thesis)
### 
### Author: Sebastian Meyer
### $Date: 2010-04-23 10:49:46 +0200 (Fri, 23 Apr 2010) $
###
### PARAMS:
### polyregion: a polygon of class "gpc.poly" (gpclib::) or "owin" (spatstat::)
### f:  two-dimensional integrand. The function must take a coordinate matrix as
###     its first argument.
### ...: further arguments passed to f.
################################################################################


### Wrapper function for the various cubature methods
# ...: arguments of f or of the specific methods

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



### Alternative 1: Two-dimensional midpoint rule (using spatstat image)
# polyregion can be anything coercible to "owin" by as.owin
# eps: width and height of the pixels (squares)
# dimyx: number of subdivisions in each dimension

polyCub.midpoint <- function (polyregion, f, ..., eps = NULL, dimyx = NULL, plot = FALSE)
{
    # as.im needs seperate x and y arguments
    fxy <- function (x, y, ...) f(cbind(x,y), ...)

    # calculate pixel values of fxy
# hoehle - 10 Apr 2011 - problem in new spatstat, if eps too large the
# try does not work anymore and the function interrupts. Hence, a crude
# bug fix is implemented forcing an upper limit on eps, until Sebastian
# finds time to study the code more carefully.
    #eps <- max(exp(1),min(exp(5),eps)) #force upper & limit limit on eps (this is a hack, but as.mask appears to have a problem that makes try/catch useless.
# sebastian - 23 Apr 2012 - these hidden limits are dangerous
# and i cannot observe any problem with spatstat here.
# so we just let the function do its work... (fingers crossed)
    IM <- tryCatch(
          spatstat::as.im.function(X = fxy, W = polyregion, ...,
                                   eps = eps, dimyx = dimyx),
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
    spatstat::plot.im(IM, axes = TRUE, col=grey(31:4/35), main="")
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



### Alternative 2: Product Gauss cubature of two-dimensional functions
### over simple polygons proposed by Sommariva & Vianello (2007)
# polyregion can be anything coercible to "gpc.poly", BUT:
#       If supplied as a "gpc.poly" be aware of the following restriction:
#       Vertices are assumed to be ordered according to the "sp" convention,
#       i.e. _clockwise_ for normal boundaries and _anticlockwise_ for holes.
#       However, in contrast to "sp", the first vertex should NOT be repeated!
#       The coerce-methods from "Polygon" or "owin" to "gpc.poly" defined in
#       package surveillance follow this convention for "gpc.poly".
# N: number of nodes of 1-dimensional Gauss-Legendre rule (see gaussCub)
# a: base-line at x = a (see gaussCub)
#    a = 0 seems to be reasonable if f has its maximum value at (0,0), e.g. for
#    the bivariate normal density with zero mean

polyCub.SV <- function (polyregion, f, ..., N, a = 0, plot = FALSE)
{
    if (!require("statmod")) {
        stop("package ", sQuote("statmod"), " is needed for Gaussian cubature")
    }
    polyregion <- as(polyregion, "gpc.poly")
    polys <- polyregion@pts
    if (plot) {
    	plot(polyregion, poly.args=alist(lwd=2), ann=FALSE)
        for (i in seq_along(polys)) polys[[i]]$ID <- i
    }
    
    respolys <- sapply(polys, function(poly) {
        # gaussCub function assumes anticlockwise order
        # (clockwise order leads to inverse integral value needed for holes)
        x <- rev(poly$x)
        y <- rev(poly$y)
        nw <- gaussCub(c(x,x[1]), c(y,y[1]), N = N, a = a)
        if (plot) points(nw$nodes, cex=0.6, pch = poly$ID) #, col=1+(nw$weights<=0)
        fvals <- f(nw$nodes, ...)
        cubature_val <- sum(nw$weights * fvals)
# if (!isTRUE(all.equal(0, cubature_val))) {
# if ((1 - 2 * as.numeric(poly$hole)) * sign(cubature_val) == -1)
# warning("wrong sign if positive integral")
# }
        cubature_val
    })
    int <- sum(respolys)
    int
}



### Function to calculate nodes and weights of the product Gauss cubature
### Code is based on the MATLAB implementation of Sommariva & Vianello (2007)
# Note: The efficiency rotation also proposed by Sommariva & Vianello (2007)
# is not implemented since it does in general only apply to convex polygons.
# Potential improvement: The efficient implementation of this function in C
# seems possible and would increase the speed of the cubature
# Parameters:
# x_bd, y_bd: coordinates of the polygon's vertices in _anticlockwise_ order
#             (otherwise the result of the cubature will have a negative sign)
#             with _repeated first vertex_ at the end
# N: degree of the one-dimensional Gauss-Legendre quadrature rule
#    (see statmod::gauss.quad, on which this function depends)
# a: base-line at x = a (see the referenced paper for an explication).
#    If NULL (the default), the midpoint of the x-range is chosen.

gaussCub <- function (x_bd, y_bd, N = 10, a = NULL)
{
    # NUMBER OF SIDES OF THE POLYGON.
    L <- length(x_bd) - 1L
    
    # base-line at x=a
    if (is.null(a)) {
        xrange <- range(x_bd)
        a <- (xrange[1] + xrange[2]) / 2
    }
    
    
    # %-------------------------------------------------------------------------
    # % COMPUTE NODES AND WEIGHTS OF 1D GAUSS-LEGENDRE RULE.
    # %-------------------------------------------------------------------------
    
    # % DEGREE "N" (as requested) (ORDER GAUSS PRIMITIVE)
    nw_N <- statmod::gauss.quad(n = N, kind = "legendre")

    # # % DEGREE "M" = N+1 (ORDER GAUSS INTEGRATION)
    M <- N + 1L
    nw_M <- statmod::gauss.quad(n = M, kind = "legendre")
    
    
    # %-------------------------------------------------------------------------
    # % COMPUTE 2D NODES (nodes_x,nodes_y) AND WEIGHTS "weights".
    # %-------------------------------------------------------------------------
    
    maxRows <- L*M
    nodes_x <- rep.int(0, maxRows*N)
    dim(nodes_x) <- c(maxRows, N)
    nodes_y <- nodes_x
    weights <- numeric(maxRows)
    K <- 0L
    
    for (side in seq_len(L))
    {
        x1 <- x_bd[side];  x2 <- x_bd[side+1L]
        y1 <- y_bd[side];  y2 <- y_bd[side+1L]
        
        if ((x1 == a && x2 == a) || (y2 == y1)) {
            # side lies on base-line or is orthogonal to it
            next
        }

        if (x2 == x1) { # side is parallel to base-line => degree N
            n_loc <- nw_N$nodes
            w_loc <- nw_N$weights
        } else { # degree M=N+1
            n_loc <- nw_M$nodes
            w_loc <- nw_M$weights
        }
        degree_loc <- length(n_loc)   # = N or M
        
        half_pt_x <- (x1+x2)/2;  half_length_x <- (x2-x1)/2;
        half_pt_y <- (y1+y2)/2;  half_length_y <- (y2-y1)/2;
        
        #% GAUSSIAN POINTS ON THE SIDE.
        x_gauss_side <- half_pt_x + half_length_x * n_loc
        y_gauss_side <- half_pt_y + half_length_y * n_loc
        
#         scaling_fact <- x_gauss_side / 2
#         scaling_fact_plus <- (x_gauss_side + a) / 2
        scaling_fact_minus <- (x_gauss_side - a) / 2
            
        local_weights <- (half_length_y * scaling_fact_minus) * w_loc
        
#         term_1 <- matrix(scaling_fact_plus, nrow = degree_loc, ncol = N)
        term_2 <- matrix(scaling_fact_minus, nrow = degree_loc, ncol = N)
#         term <- matrix(scaling_fact, nrow = degree_loc, ncol = N)
        rep_n_Np1 <- matrix(nw_N$nodes + 1, nrow = degree_loc, ncol = N, byrow = TRUE)
        
        #% x, y ARE STORED IN MATRICES. A COUPLE WITH THE SAME INDEX IS A POINT,
        #% i.e. "P_i=(x(k),y(k))" FOR SOME "k".
        # x = (term+a/2) + (term-a/2) * rep_n_N = (1+rep_n_N) * term + (1-rep_n_N) * a/2
        # x = (term_2+a) + term_2 * rep_n_N = (1+rep_n_N) * term_2 + a
        x <- rep_n_Np1 * term_2 + a
#         y = matrix(y_gauss_side, nrow = degree_loc, ncol = N)
        
        # add nodes and weights for this side of the polygon
        rowidx <- K + seq_len(degree_loc)
        nodes_x[rowidx,] <- x
        nodes_y[rowidx,] <- y_gauss_side   # fills (sub-)matrix by column
        weights[rowidx] <- local_weights
        K <- K + degree_loc
    }
    
    # only the first K entries of 'weights' and rows of 'nodes_x' and 'nodes_y'
    # have been filled, the remainder till 'maxRows', contains the zeros
    # from initialisation.
    if (K < maxRows) {
        seqK <- seq_len(K)
        weights <- weights[seqK]
        nodes_x <- nodes_x[seqK,]
        nodes_y <- nodes_y[seqK,]
    }
    nodes <- cbind(c(nodes_x), c(nodes_y))   # K*N x 2 matrix
    weightsvec <- rep(nw_N$weights, each = K) * rep.int(weights, N)  # K*N long
#     f_xy <- f(nodes, ...)   # K*N function evaluations
#     dim(f_xy) <- c(K, N)
#     .tmp <- colSums(weights * f_xy)   # equals t(weights) %*% f_xy, which is 1 x N
#     cubature_val <- sum(.tmp * nw_N$weights)   # equals .tmp %*% nw_N$weights

    ret <- list(nodes = nodes, weights = weightsvec)
    return(ret)
}



### Alternative 3: Quasi-exact cubature specific for the bivariate normal density
### based on triangulation and formulae from Chapter 26 of the famous
### Abramowitz & Stegun handbook (cf. Section 26.9, Example 9, pp. 956f. therein)
# quite complicated because A&S formula is only for triangles where one vertex is the origin (0,0)
# For each triangle of the tristrip we have to check in which of the 6 outer regions
# of the triangle the origin (0,0) lies and adapt the signs in the formula adequately
# (AOB+BOC-AOC) or (AOB-AOC-BOC) or (AOB+AOC-BOC) or (AOC+BOC-AOB) or ...
# In contrast, the most time consuming step is the evaluation of pmvnorm in .V

# helper function, which calculates the integral of the standard bivariat normal
# over a triangle bounded by y=0, y=ax, x=h (cf. formula 26.3.23)
.V <- function(h,k) {
    a <- k/h
    rho <- -a/sqrt(1+a^2)
    # V = 0.25 + L(h,0,rho) - L(0,0,rho) - Q(h) / 2
    # L(0,0,rho) = 0.25 + asin(rho) / (2*pi)
    # V = L(h,0,rho) - asin(rho)/(2*pi) - Q(h) / 2
    Lh0rho <- mvtnorm::pmvnorm(
        lower = c(h,0), upper = c(Inf,Inf), mean = c(0,0), corr = matrix(c(1,rho,rho,1),2,2)
    )
    Qh <- stats::pnorm(h, mean = 0, sd = 1, lower.tail = FALSE)
    return(Lh0rho - asin(rho)/2/pi - Qh/2)
}

# helper function, which calculates the integral of the standard bivariat normal over a triangle A0B
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

# checks if point1 and point2 lie on the same side of a line through linepoint1 and linepoint2
# see http://www.gamedev.net/community/forums/topic.asp?topic_id=457450
.pointsOnSameSide <- function (linepoint1, linepoint2, point1, point2 = c(0,0))
{
    n <- c(-1,1) * rev(linepoint2-linepoint1)   # normal vector
    S <- dotprod(point1-linepoint1,n) * dotprod(point2-linepoint1,n)
    return(S > 0)
}

# helper function, which calculates the integral of the standard bivariat normal over a triangle ABC
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

polyCub.exact.Gauss <- function (polyregion, mean = c(0,0), Sigma = diag(2), plot = FALSE)
{
	if (!require("mvtnorm")) {
        stop("package ", sQuote("mvtnorm"), " is needed for the exact Gaussian method")
    }
    polyregion <- as(polyregion, "gpc.poly")
    
    # coordinate transformation so that the standard bivariat normal density
    # can be used in integrations (cf. formula 26.3.22)
    rho <- stats::cov2cor(Sigma)[1,2]
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

