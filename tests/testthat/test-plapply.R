context("plapply()")

test_that("plapply() results are reproducible", {
    res1 <- plapply(c(1, 1), rnorm, .parallel = 2, .seed = 1, .verbose = FALSE)
    res2 <- plapply(c(1, 1), rnorm, .parallel = 2, .seed = 1, .verbose = FALSE)
    expect_identical(res1, res2)
})
