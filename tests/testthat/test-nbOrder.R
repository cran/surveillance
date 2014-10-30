context("Neighbourhood order")

## generate random adjancency matrix
set.seed(3)
n <- 5
adjmat <- matrix(0, n, n, dimnames=list(letters[1:n],letters[1:n]))
adjmat[lower.tri(adjmat)] <- sample(0:1, n*(n-1)/2, replace=TRUE)
adjmat <- adjmat + t(adjmat)

## validated matrix of neighbourhood orders
nbmat <- structure(
    c(0L, 2L, 1L, 3L, 2L,
      2L, 0L, 1L, 1L, 2L,
      1L, 1L, 0L, 2L, 1L,
      3L, 1L, 2L, 0L, 1L,
      2L, 2L, 1L, 1L, 0L),
    .Dim = c(5L, 5L),
    .Dimnames = rep.int(list(c("a", "b", "c", "d", "e")), 2L)
    )

if (requireNamespace("spdep")) {
    test_that("nbOrder() returns the validated matrix", {
        expect_that(nbOrder(adjmat, maxlag=Inf),
                    is_identical_to(nbmat))
    })
}

test_that("zetaweights(.,maxlag=1,normalize=FALSE) is inverse of nbOrder", {
    expect_that(zetaweights(nbmat, maxlag=1, normalize=FALSE),
                is_identical_to(adjmat))
})
