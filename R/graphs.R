################################################################################
### Functions concerning graphs: neighbourhood order, adjacency matrix
###
### Copyright (C) 2009-2013,2017,2023 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################


### Determine the matrix of neighbourhood orders
### given the binary matrix of first-order neighbours.

nbOrder <- function (neighbourhood, maxlag = Inf)
{
    stopifnot(isScalar(maxlag), maxlag > 0)
    checkNeighbourhood(neighbourhood)
    neighbourhood <- neighbourhood == 1           # convert to binary matrix
    nregions <- nrow(neighbourhood)
    maxlag <- as.integer(min(maxlag, nregions-1)) # upper bound of nb order

    if (maxlag == 1L) {
        storage.mode(neighbourhood) <- "integer"
        return(neighbourhood)
    }

    ## list of indexes of first-order neighbours by region
    ##first <- apply(neighbourhood, 1L, which, simplify = FALSE) # R >= 4.1.0
    first <- lapply(seq_len(nregions), function (i) which(neighbourhood[i,]))

    ## Side note: fast method to determine neighbours _up to_ specific order:
    ## crossprod(neighbourhood) > 0  # up to second order neighbours (+set diag to 0)
    ## (neighbourhood %*% neighbourhood %*% neighbourhood) > 0  # up to order 3
    ## and so on...

    ## now find recursive neighbours for each region
    nbmat <- `diag<-`(neighbourhood, NA)  # skip self in which() below
    for (i in seq_len(nregions)) {  # slightly faster than [l]apply() variants
        nblags <- as.integer(nbmat[i,])
        lag <- 1L
        while (lag < maxlag) {
            nbs <- which(nblags == lag)
            nbs2 <- unlist(first[nbs])
            new <- intersect(nbs2, which(nblags == 0))
            if (length(new)) {
                lag <- lag + 1L
                nblags[new] <- lag
            } else break
        }
        nbmat[i,] <- nblags
    }
    diag(nbmat) <- 0L

    ## Done
    nbmat
}


### Derive adjacency structure from a SpatialPolygons object
### Working horse: spdep::poly2nb

poly2adjmat <- function (SpP, ..., zero.policy = TRUE)
{
    if (!requireNamespace("spdep"))
        stop("package ", sQuote("spdep"),
             " is required to derive adjacencies from SpatialPolygons")
    nb <- spdep::poly2nb(SpP, ...)
    adjmat <- spdep::nb2mat(nb, style="B", zero.policy=zero.policy)
    attr(adjmat, "call") <- NULL
    colnames(adjmat) <- rownames(adjmat)
    adjmat
}
