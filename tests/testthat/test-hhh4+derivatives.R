### Fixed effects hhh4() model fit and involved analytical derivatives

data("measlesWeserEms")
measlesModel <- list(
    end = list(f = addSeason2formula(~1 + t, S=1, period=52),
               offset = population(measlesWeserEms)),
    ar = list(f = ~1),
    ne = list(f = ~1 + log(pop),
        weights = W_powerlaw(maxlag = 5, normalize = TRUE)),
    family = "NegBin1", data = list(pop = population(measlesWeserEms))
)

measlesFit <- hhh4(stsObj = measlesWeserEms, control = measlesModel)

test_that("estimates and standard errors are reproducible", {
    ## dput(coef(measlesFit, se = TRUE))
    orig <- structure(
        c(-0.499636482022272, 0.551345030080107, 0.96093157194767,
          -0.153585641356373, 0.00333284018297979, 1.01500011496702,
          -0.588738943313705, 5.52782609236691, 1.81915612994789,
          0.121781347106564, 1.27401298230559, 0.453889365025671,
          0.281013375484401, 0.00459840327748742, 0.210642721317572,
          0.191921649336323, 1.87984346848385, 0.265016986696184),
        .Dim = c(9L, 2L),
        .Dimnames = list(c("ar.1", "ne.1", "ne.log(pop)", "end.1",
            "end.t", "end.sin(2 * pi * t/52)", "end.cos(2 * pi * t/52)",
            "neweights.d", "overdisp"), c("Estimate", "Std. Error"))
    )
    expect_equal(coef(measlesFit, se = TRUE), orig,
                 tolerance = 1e-6) # increased for Solaris Sparc
    ## tolerance determined empirically by an R build with --disable-long-double
})

test_that("neighbourhood weights array yields the same results", {
    What <- getNEweights(measlesFit)
    ## put that in an array for time-varying weights in hhh4
    ## (they are not actually varying here)
    Warray <- array(What,
                    dim = c(dim(What),nrow(measlesWeserEms)),
                    dimnames = c(dimnames(What), list(NULL)))
    measlesFit_Warray <- update(measlesFit, ne = list(weights = Warray),
                                use.estimates = FALSE)
    ## NOTE: variance estimates are different because of fixed powerlaw
    expect_equal(measlesFit_Warray, measlesFit,
                 ignore = c("control", "coefficients", "se", "cov", "dim"))
    expect_equal(coef(measlesFit_Warray),
                 coef(measlesFit)[names(coef(measlesFit_Warray))],
                 tolerance = 1e-6)  # triggered by 64-bit win-builder
})

test_that("score vector and Fisher info agree with numerical approximations", if (requireNamespace("numDeriv")) {
    test <- function (neweights) {
        Wname <- deparse(substitute(neweights))
        measlesModel$ne$weights <- neweights
        capture.output( # hide reports as we use a different tolerance
        pencomp <- hhh4(measlesWeserEms, measlesModel,
                        check.analyticals = "numDeriv")$pen
        )
        expect_equal(pencomp$score$analytic, pencomp$score$numeric,
                     tolerance = .Machine$double.eps^0.5, info = Wname)
        expect_equal(pencomp$fisher$analytic, pencomp$fisher$numeric,
                     tolerance = .Machine$double.eps^0.25, info = Wname)
    }
    test(W_powerlaw(maxlag = 5, normalize = FALSE, log = FALSE))
    ## normalized PL with maxlag < max(nbmat) failed in surveillance < 1.9.0:
    test(W_powerlaw(maxlag = 3, normalize = TRUE, log = TRUE))
    ## check unconstrained weights
    test(W_np(maxlag = 5, truncate = TRUE, normalize = FALSE))
    test(W_np(maxlag = 3, truncate = FALSE, normalize = TRUE))
    ## test two-component formulations (AR within NE)
    measlesModel$ar <- list(f = ~ -1)
    test(W_powerlaw(maxlag = 3, normalize = TRUE, log = TRUE, from0 = TRUE))
    test(W_np(maxlag = 1, truncate = FALSE, normalize = FALSE, from0 = TRUE))
    test(W_np(maxlag = 3, truncate = TRUE, normalize = TRUE, from0 = TRUE))
})

test_that("alternative optimizers give equivalent results", {
    ctrl_nlm <- list(method = "nlm", check.analyticals = TRUE)
    expect_equal(
        update(measlesFit, optimizer = list(regression = ctrl_nlm),
               use.estimates = FALSE),
        measlesFit, tolerance = 1e-5, ignore = "control"
    )
    ctrl_BFGS <- list(method = "BFGS", reltol = 1e-12)
    expect_equal(
        update(measlesFit, optimizer = list(regression = ctrl_BFGS),
               use.estimates = FALSE),
        measlesFit, tolerance = 1e-5, ignore = "control"
    )
})

test_that("automatic and manual normalization are equivalent", {
    ## check for equivalent functions
    for (type in c("powerlaw", "np")) {
        W_type <- get(paste0("W_", type), mode = "function")
        w0 <- W_type(maxlag = 3, normalize = TRUE)
        w1 <- surveillance:::scaleNEweights.list(
            W_type(maxlag = 3, normalize = FALSE),
            normalize = TRUE)
        pars <- w0$initial
        nbmat <- neighbourhood(measlesWeserEms)
        expect_equal(w1$w(pars, nbmat), w0$w(pars, nbmat))
        ## for the power law, dw and d2w are length 1 lists in w1 but not in w0
        unlistIfPL <- if (type == "powerlaw") function (x) x[[1L]] else identity
        expect_equal(unlistIfPL(w1$dw(pars, nbmat)), w0$dw(pars, nbmat))
        expect_equal(unlistIfPL(w1$d2w(pars, nbmat)), w0$d2w(pars, nbmat))
        ## microbenchmark::microbenchmark(w1$d2w(pars, nbmat), w0$d2w(pars, nbmat))
        ## -> type-specific implementations of normalized derivatives are faster
    }
    ## check for equivalent fits (rather redundant)
    measlesFit2 <- hhh4(
        stsObj = measlesWeserEms,
        control = modifyList(measlesModel, list(
            ne = list(
                weights = W_powerlaw(maxlag = 5, normalize = FALSE),
                normalize = TRUE # -> use scaleNEweights.list()
                )))
        )
    expect_equal(measlesFit, measlesFit2, ignore = "control",
                 tolerance = 1e-6) # increased to pass on 32-bit Windows
})

test_that("unnamed plot() 'type' argument is not passed down", {
    expect_length(plot(measlesFit, "neweights", plotter = list), 1)
    ## failed in surveillance <= 1.24.1; actual use case:
    ## R> plot(measlesFit, "neweights", exclude = NULL)
    ## Error in get(as.character(FUN), mode = "function", envir = envir) :
    ##   object 'neweights' of mode 'function' was not found
})

measlesWeserEms2 <- measlesWeserEms
neighbourhood(measlesWeserEms2) <- neighbourhood(measlesWeserEms2) + 1L

test_that("W_powerlaw(..., from0 = TRUE) equals manual approach", {
    measlesModel2 <- modifyList(measlesModel, list(
        ar = list(f = ~ -1),
        ne = list(weights = W_powerlaw(maxlag = 5, from0 = TRUE))
        ))
    measlesFit2 <- hhh4(measlesWeserEms, measlesModel2)
    ## manual approach
    measlesModel2_manual <- modifyList(measlesModel2, list(
        ne = list(weights = W_powerlaw(maxlag = 5 + 1))
        ))
    measlesFit2_manual <- hhh4(measlesWeserEms2, measlesModel2_manual)
    expect_equal(measlesFit2, measlesFit2_manual,
                 ignore = c("control", "stsObj"))
})

test_that("W_np(..., from0 = TRUE) equals manual approach", {
    measlesModel2 <- modifyList(measlesModel, list(
        ar = list(f = ~ -1),
        ne = list(weights = W_np(maxlag = 2, from0 = TRUE))
        ))
    measlesFit2 <- hhh4(measlesWeserEms, measlesModel2)
    ## manual approach
    measlesModel2_manual <- modifyList(measlesModel2, list(
        ne = list(weights = W_np(maxlag = 2 + 1))
        ))
    measlesFit2_manual <- hhh4(measlesWeserEms2, measlesModel2_manual)
    expect_equal(measlesFit2, measlesFit2_manual,
                 ignore = c("control", "stsObj"))
})
