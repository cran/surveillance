context("hhh4() with epidemic offsets")

data("measlesWeserEms")
measles2 <- measlesWeserEms[,c("03457","03454")]

## AR model
fit1 <- hhh4(measles2, list(ar = list(f = ~1)))
##plot(fit1, units=NULL)

## use estimated exp(lambda) as offset -> new lambda should be 0, equal fit
o1 <- exp(fit1$coefficients[["ar.1"]])
fit1o <- hhh4(measles2, list(
    ar = list(f = ~1, offset = matrix(o1, nrow(measles2), ncol(measles2)))
    ))
test_that("model with AR offset is fitted correctly", {
    expect_equal(fit1o$coefficients[["ar.1"]], 0)
    expect_equal(fitted(fit1o), fitted(fit1))
})

## same test with an AR+NE model
fit2 <- hhh4(measles2, list(ar = list(f = ~1), ne = list(f = ~1)))
##plot(fit2, units=NULL)
o2_ar <- exp(fit2$coefficients[["ar.1"]])
o2_ne <- exp(fit2$coefficients[["ne.1"]])
fit2o <- hhh4(measles2, list(
    ar = list(f = ~1, offset = matrix(o2_ar, nrow(measles2), ncol(measles2))),
    ne = list(f = ~1, offset = matrix(o2_ne, nrow(measles2), ncol(measles2)))
    ))
test_that("model with AR+NE offsets is fitted correctly", {
    expect_equal(fit2o$coefficients[["ar.1"]], 0, scale = 1)  # use abs. diff
    expect_equal(fit2o$coefficients[["ne.1"]], 0, scale = 1,
                 tolerance = 1e-6)  # for ATLAS/MKL/OpenBLAS
    expect_equal(fitted(fit2o), fitted(fit2))
})

## createLambda() and thus maxEV was wrong in surveillance <= 1.16.1
test_that("Lambda matrix incorporates epidemic offsets", {
    expect_equal(getMaxEV(fit1o)[1], getMaxEV(fit1)[1])
    expect_equal(getMaxEV(fit2o)[1], getMaxEV(fit2)[1])
})

## simulate.hhh4() was wrong in surveillance <= 1.16.1
test_that("simulation accounts for epidemic offsets", {
    ## check the relative difference in the total number of cases
    obs <- fitted(fit2o)
    sim <- simulate(fit2o, seed = 1, y.start = observed(measles2)[1,],
                    subset = fit2o$control$subset, simplify = TRUE)
    expect_true(abs(sum(sim)/sum(obs)-1) < 0.5)
})
