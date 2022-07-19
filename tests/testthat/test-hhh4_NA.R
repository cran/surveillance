### fitting hhh4() models to time series with missing values

data("influMen")
fluMen <- disProg2sts(influMen)
## set some observations to NA
set.seed(3)
is.na(fluMen@observed) <- sample(length(fluMen@observed), 100)

## compare endemic-only model against NegBin-GLM
form <- addSeason2formula(f = ~ -1 + fe(1, which = c(TRUE, TRUE)), S = c(3, 1))
fitHHH <- hhh4(fluMen,
               list(end = list(f=form), family = "NegBin1", subset = 1:nrow(fluMen)))
fitGLM <- MASS::glm.nb(
    formula = observed ~ -1 + unit + sin(2*pi*t/52):unit + cos(2*pi*t/52):unit +
        I(sin(4*pi*t/52)*unitI) + I(cos(4*pi*t/52)*unitI) +
        I(sin(6*pi*t/52)*unitI) + I(cos(6*pi*t/52)*unitI),
    data = transform(tidy.sts(fluMen), t = epoch - 1, unitI = unit == "influenza"))

expect_equal(logLik(fitHHH), logLik(fitGLM))
expect_equal(fitted(fitHHH)[!terms(fitHHH)$isNA], unname(fitted(fitGLM)))
expect_equivalent(coef(fitHHH)[["overdisp"]], 1/fitGLM$theta)
idxhhh <- c(1:2, 7:10, 3:6)
expect_equivalent(head(fitHHH$coefficients, -1), fitGLM$coefficients[idxhhh])
expect_equivalent(head(fitHHH$se, -1), summary(fitGLM)$coefficients[idxhhh, 2],
                  tolerance = 0.01)

### compare AR-only model against NegBin-GLM
## meningococcal counts are strictly positive so plain AR works
men <- fluMen[,"meningococcus"]
fitHHH_AR <- hhh4(men, list(end = list(f = ~-1), family = "NegBin1",
                            ar = list(f = addSeason2formula(~1))))
fitGLM_AR <- MASS::glm.nb(
    formula = addSeason2formula(observed ~ 1 + offset(log(Ylag))),
    data = transform(tidy.sts(men), t = epoch - 1, Ylag = c(NA, head(observed, -1))))

expect_equal(logLik(fitHHH_AR), logLik(fitGLM_AR))
expect_equal(fitted(fitHHH_AR)[!terms(fitHHH_AR)$isNA[fitHHH_AR$control$subset,,drop=FALSE]],
             unname(fitted(fitGLM_AR)))
expect_equivalent(coef(fitHHH_AR)[["overdisp"]], 1/fitGLM_AR$theta)
expect_equivalent(head(fitHHH_AR$coefficients, -1), fitGLM_AR$coefficients)
expect_equivalent(head(fitHHH_AR$se, -1), summary(fitGLM_AR)$coefficients[, 2],
                  tolerance = 0.05)

### compare NE-only model against NegBin-GLM (where NE is actually AR as above)
expect_warning(
fitHHH_NE <- hhh4(men, list(end = list(f = ~-1), family = "NegBin1",
                            ne = list(f = addSeason2formula(~1), weights = diag(1))))
, "requires a multivariate")
expect_equivalent(fitHHH_AR, fitHHH_NE, ignore = c("control", "lags"))
expect_equal(meanHHH(fitHHH_AR$coefficients, terms(fitHHH_AR))$epi.own,
             meanHHH(fitHHH_NE$coefficients, terms(fitHHH_NE))$epi.neighbours)
if (dev.capabilities("capture")[[1L]]) { # e.g. not in tinytest as that uses pdf
plot(fitHHH_AR, legend = FALSE, col = c(8,8,8)); plotARfit <- dev.capture()
plot(fitHHH_NE, legend = FALSE, col = c(8,8,8)); plotNEfit <- dev.capture()
expect_identical(plotARfit, plotNEfit)
}
