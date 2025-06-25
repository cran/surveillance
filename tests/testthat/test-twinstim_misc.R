## miscellaneous regression tests for twinstim()

data("imdepi")
load(system.file("shapes", "districtsD.RData", package = "surveillance"))


## let some districts have no population
imdepi0 <- update(imdepi,
    stgrid = within(imdepi$stgrid,
                    popdensity[startsWith(as.character(tile), "01")] <- 0))
stopifnot(# verify update-method
    identical(names(imdepi0$stgrid), names(imdepi$stgrid)),
    identical(names(imdepi0$events), names(imdepi$events)),
    with(imdepi0$events@data,
         all.equal(startsWith(as.character(tile), "01"), popdensity == 0)),
    identical(marks(imdepi0), marks(imdepi))
)

## automatic start value is robust against -Inf offset
tools::assertWarning(# infinite logLik
fit0 <- twinstim(endemic = ~offset(log(popdensity)) + I(start/365),
                 data = imdepi0, model = TRUE,
                 optim.args = list(fixed = TRUE), verbose = FALSE)
)
## beta0 was initialized at Inf in surveillance <= 1.22.1
stopifnot(is.finite(coef(fit0)), is.character(fit0$converged),
          is.infinite(logLik(fit0))) # because of events in 0-pop tiles

## endemic intensity is 0 in unpopulated districts
hGrid <- intensity.twinstim(fit0, "space", tiles = districtsD)$hGrid
districts_h0 <- names(which(hGrid == 0))
stopifnot(length(districts_h0) > 0, startsWith(districts_h0, "01"))


## intensityplot works for an endemic-only model
fit_end <- twinstim(endemic = ~1, data = imdepi, model = TRUE,
                    optim.args = list(fixed = TRUE), verbose = FALSE)
intensityplot(fit_end, "total", "space", tiles = districtsD) -> .plotobj
## produced an error in surveillance <= 1.22.1:
##   unable to find an inherited method for function 'coordinates' for signature '"NULL"'
