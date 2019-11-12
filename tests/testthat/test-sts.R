context("S4 class definition of \"sts\" and its extensions")

test_that("\"sts\" prototype is a valid object",
          expect_true(validObject(new("sts"))))

mysts <- sts(1:10, frequency = 4, start = c(1959, 2))

test_that("conversion from \"ts\" to \"sts\" works as expected", {
    myts <- ts(1:10, frequency = 4, start = c(1959, 2))
    expect_identical(as(myts, "sts"), mysts)
    ## this failed in surveillance 1.11.0 due to a wrong "start" calculation
})

test_that("if missing(observed), initialize-method copies slots", {
    mysts_updated <- initialize(mysts, epoch = 2:11)
    expect_identical(mysts_updated@epoch, 2:11)
    mysts_updated@epoch <- mysts@epoch
    expect_identical(mysts_updated, mysts)
    ## construct stsBP from existing "sts" object
    mystsBP <- new("stsBP", mysts,
                   ci = array(NA_real_, c(10,1,2)),
                   lambda = array(NA_real_, c(10,1,1)))
    expect_identical(as(mystsBP, "sts"), mysts)
})

test_that("different initializations of \"stsBP\" work as expected", {
    mystsBP <- new("stsBP", observed = 1:10, freq = 4, start = c(1959, 2),
                   ci = array(NA_real_, c(10,1,2)),
                   lambda = array(NA_real_, c(10,1,0)))
    expect_identical(mystsBP, as(mysts, "stsBP"))
})

test_that("different initializations of \"stsNC\" work as expected", {
    mystsNC <- new("stsNC", observed = 1:10, freq = 4, start = c(1959, 2),
                   pi = array(NA_real_, c(10,1,2)),
                   SR = array(NA_real_, c(10,0,0)))
    expect_identical(mystsNC, as(mysts, "stsNC"))
})

test_that("sts(..., population) sets the populationFrac slot", {
    ## for sts() construction, "population" is an alias for "populationFrac"
    ## (the internal slot name), introduced in the space-time JSS paper
    sts1 <- sts(cbind(1:3, 11:13), population = c(10, 20))
    sts2 <- sts(cbind(1:3, 11:13), populationFrac = c(10, 20))
    expect_identical(sts1, sts2)
})

test_that("\"sts\" conversion to a (tidy) data frame works consistently", {
    ## univariate sts
    mystsdata <- as.data.frame(mysts, as.Date = FALSE)
    expect_identical(tidy.sts(mysts)[names(mystsdata)], mystsdata)
    ## multivariate sts
    data("momo")
    momo3tidy_uv <- tidy.sts(momo[,3])
    momo3tidy_mv <- subset(tidy.sts(momo), unit == levels(unit)[3])
    momo3tidy_mv$unit <- momo3tidy_mv$unit[drop=TRUE]
    row.names(momo3tidy_mv) <- NULL
    expect_identical(momo3tidy_uv, momo3tidy_mv)
})

test_that("we can subset epochs of an \"sts\" object", {
    expect_identical(mysts[TRUE,TRUE], mysts)
    expect_identical(mysts[2,]@start, c(1959, 3))
    ## negative and 0 indices produced wrong "start" in surveillance <= 1.16.2
    expect_identical(mysts[-1,], mysts[2:10,])
    expect_identical(mysts[0,]@start, mysts@start)
})

test_that("colnames need to be identical (only for multivariate data)", {
    slots_dn <- c("observed", "state", "alarm", "upperbound", "populationFrac")
    ## ignore colnames mismatch for univariate time series
    sts_args_1 <- lapply(setNames(nm = slots_dn), function (slot)
        matrix(0, 1, 1, dimnames = list(NULL, slot)))
    sts_args_1$neighbourhood <- matrix(0, 1, 1, dimnames = list("a", "a"))
    expect_silent(do.call(sts, sts_args_1))
    ## multivariate time series with inconsistent column order are invalid
    sts_args_2 <- list(
        observed = matrix(0, 1, 2, dimnames = list(NULL, c("r1", "r2")))
    )
    sts_args_2[slots_dn[-1]] <- list(sts_args_2$observed[,2:1,drop=FALSE])
    sts_args_2$neighbourhood <- matrix(0, 2, 2, dimnames = rep(list(c("r2", "r1")), 2))
    expect_error(do.call(sts, sts_args_2), "colnames")  # new in surveillance > 1.17.1
    ## column names can be missing for other slots
    expect_silent(do.call(sts, c(sts_args_2[1], lapply(sts_args_2[-1], unname))))
})
