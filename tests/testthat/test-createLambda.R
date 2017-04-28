context("Create next-generation matrix Lambda from a \"hhh4\" model")

data("measlesWeserEms")

## a simple endemic model
measlesFit0 <- hhh4(measlesWeserEms, list(
    end = list(f = addSeason2formula(~1), offset = population(measlesWeserEms)),
    family = "NegBin1"
))

test_that("endemic-only model has zero-valued Lambda matrix", {
    res <- getMaxEV_season(measlesFit0)
    expect_equal(res$maxEV.const, 0)
    zeromat <- matrix(0, measlesFit0$nUnit, measlesFit0$nUnit)
    expect_equal(res$Lambda.const, zeromat)
    expect_equal(createLambda(measlesFit0)(2), zeromat)
})

## + AR component
measlesFit1 <- update(measlesFit0, ar = list(f = addSeason2formula(~1)))

test_that("autoregressive model has a diagonal Lambda matrix", {
    res <- getMaxEV_season(measlesFit1)
    expect_equal(res$Lambda.const, diag(res$maxEV.const, measlesFit1$nUnit))
    expect_equal(createLambda(measlesFit1)(2),
                 diag(res$maxEV.season[2], measlesFit1$nUnit))
})

## + NE component
measlesFit2 <- update(measlesFit1,
    ne = list(f = ~1, weights = neighbourhood(measlesWeserEms) == 1)) # symmetric
measlesFit3 <- update(measlesFit2,
    ne = list(normalize = TRUE)) # asymetric

test_that("getMaxEV() and getMaxEV_season() agree", {
    expect_equal(getMaxEV_season(measlesFit2)$maxEV.season,
                 getMaxEV(measlesFit2)[seq_len(measlesWeserEms@freq)])
    expect_equal(getMaxEV_season(measlesFit3)$maxEV.season,
                 getMaxEV(measlesFit3)[seq_len(measlesWeserEms@freq)])
})

## AR within NE + unit-specific epidemic covariate
measlesFit4 <- update(measlesFit0,
    ne = list(f = ~pop, weights = (neighbourhood(measlesWeserEms)+1)^-2,
              normalize = TRUE),
    data = list(pop = population(measlesWeserEms)))

## calculate "nu + Lambda Y_{t-1}" and compare to fitted(object)
check_createLambda <- function (object)
{
    mname <- deparse(substitute(object))
    model <- terms(object)
    means <- meanHHH(object$coefficients, model, subset = seq_len(model$nTime))
    expect_equal(means$mean[model$subset,,drop=FALSE], fitted(object),
                 expected.label = paste0("fitted(", mname, ")"))
    Lambda <- createLambda(object)
    if (any(object$lags != 1, na.rm = TRUE))
        stop("check not implemented for lags != 1")
    meansByLambda <- t(vapply(
        X = object$control$subset,
        FUN = function(t) means$endemic[t,] + Lambda(t) %*% model$response[t-1,],
        FUN.VALUE = numeric(object$nUnit), USE.NAMES = FALSE))
    expect_equal(meansByLambda, unname(fitted(object)),
                 expected.label = paste0("fitted(", mname, ")"))
}

test_that("multivariate formulation using Lambda agrees with fitted values", {
    check_createLambda(measlesFit0)
    check_createLambda(measlesFit1)
    check_createLambda(measlesFit2)
    check_createLambda(measlesFit3)  # failed in surveillance < 1.13.1
    check_createLambda(measlesFit4)  # failed in surveillance < 1.13.1
})
