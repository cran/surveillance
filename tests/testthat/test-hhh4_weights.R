### Neighbourhood weights in hhh4()

observed <- cbind(c(1,2,4), c(1,2,4))

test_that("AR-only and NE-only fit agree in toy scenario", {
    counts <- sts(observed)
    m1 <- hhh4(counts, control = list(
        end = list(f = ~ -1), family = "Poisson",
        ar = list(f = ~1)))
    expect_equivalent(coef(m1, idx2Exp=TRUE), 2)
    ## same fit via NE (because units have identical counts)
    m2 <- hhh4(counts, control = list(
        end = list(f = ~ -1), family = "Poisson",
        ne = list(f = ~1, weights = matrix(c(0,1,1,0), 2, 2))))
    m1$control <- m2$control <- m1$lags <- m2$lags <- NULL
    expect_equivalent(m1, m2)
})

test_that("time-varying NE weights align with time index of mu", {
    W <- matrix(c(0,1,1,0), 2, 2)
    Wt <- array(c(W, W, 0*W), c(dim(W), 3))  # w_jit = 0 for t=3
    off <- surveillance:::weightedSumNE(observed, Wt, lag = 1)
    expect_true(all(is.na(off[1L,])))
    expect_identical(off[3L,], c(0, 0))  # NE sum is zero at t=3
    ## failed in surveillance <= 1.18.0, where w_ji(t-1) * y_j(t-1)
    ## was calculated, whereas w_jit * y_j(t-1) was used for simulation,
    ## the latter being the desired behaviour (same time index as covariates)
})
