library("surveillance")

## library("testthat")
## context("cross-checking neighbourhood order functions")

## generate random adjancency matrix
set.seed(1)
n <- 6
adjmat <- matrix(0, n, n, dimnames=list(letters[1:n],letters[1:n]))
adjmat[lower.tri(adjmat)] <- sample(0:1, n*(n-1)/2, replace=TRUE)
adjmat <- adjmat + t(adjmat)

## test_that("zetaweights(.,maxlag=1,normalize=FALSE) is inverse of nbOrder", {
    nbmat <- nbOrder(adjmat, maxlag=Inf)
    adjmat2 <- zetaweights(nbmat, maxlag=1, normalize=FALSE)
    ## expect_that(adjmat2, is_identical_to(adjmat))
## })

stopifnot(identical(adjmat, adjmat2))
