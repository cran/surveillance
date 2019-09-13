context("Likelihood and score function of twinstim()")
## Note: derivatives of interaction functions are tested in separate files
##       we thus use the relatively fast Gaussian kernel here

data("imdepi")
model <- twinstim(
    endemic = addSeason2formula(~offset(log(popdensity)),
        S = 1, period = 365, timevar = "start"),
    epidemic = ~type,
    siaf = siaf.gaussian(),
    tiaf = tiaf.step(2),
    data = imdepi,
    optim.args = NULL, verbose = FALSE
)
theta <- c("h.(Intercept)" = -20,
           "h.sin(2 * pi * start/365)" = 0.2, "h.cos(2 * pi * start/365)" = 0.3,
           "e.(Intercept)" = -10, "e.typeC" = -0.9,
           "e.siaf.1" = 2, "e.tiaf.1" = -1)

test_that("likelihood is still the same", {
    expect_equal(model$ll(theta), -9579.65468598488)
})

test_that("score vector agrees with numerical approximation", {
    numsc <- if (surveillance.options("allExamples") && requireNamespace("numDeriv")) {
        numDeriv::grad(func = model$ll, x = theta)
    } else { # for faster --as-cran tests
        c(-321.766081898055, -17.0779781937451, -37.1712258869585,
          -21.4444934196989, -5.43080160401029, -15.085241575699,
          -20.1708323190602)
    }
    expect_equal(model$sc(theta), numsc)
})

## Note: twinstim() uses an estimate of the _expected_ Fisher information,
##       which does not necessarily agree with the negative Hessian of the ll
##       (it does asymptotically **at the MLE**)
## numfi <- -numDeriv::hessian(func = model$ll, x = theta)
## anafi <- model$fi(theta)

test_that("one-parameter power law agrees with more general implementation", {
    m0 <- update.default(model, siaf = siaf.powerlaw(), tiaf = NULL, subset = time < 30)
    m1 <- update.default(m0, siaf = siaf.powerlaw1(sigma = exp(2)))
    expect_equal(m0$ll(theta), m1$ll(c(head(theta, -2), -1)))
    expect_equal(m0$sc(theta)[-6], m1$sc(c(head(theta, -2), -1)))
})


### now check with identity link for the epidemic predictor

model2 <- update.default(model, siaf = NULL, tiaf = NULL, epidemic = ~1, epilink = "log")
model2i <- update.default(model2, epilink = "identity")
theta2 <- theta2i <- theta[1:4]
theta2i["e.(Intercept)"] <- exp(theta2["e.(Intercept)"])

test_that("likelihoods with log-link and identity link are the same", {
    expect_equal(model2i$ll(theta2i), model2$ll(theta2))
})

test_that("identity link score vector agrees with numerical approximation", {
    numsc <- if (surveillance.options("allExamples") && requireNamespace("numDeriv")) {
        numDeriv::grad(func = model2i$ll, x = theta2i)
    } else { # for faster --as-cran tests
        c(-679.706275919901, -91.0659401491325, -114.082117122738,
          -1532144485.45524)
    }
    expect_equal(model2i$sc(theta2i), numsc)
})
