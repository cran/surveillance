## moved from example(simulate.twinstim)

data("imdepi", "imdepifit")
load(system.file("shapes", "districtsD.RData", package="surveillance"))

## simulate 2 realizations over a short period (for speed),
## with modified coefficients,
## considering original events before t0=31 as prehistory
expect_stdout(
    mysims <- simulate(imdepifit, nsim=2, seed=1, data=imdepi, tiles=districtsD,
                       newcoef=c("e.typeC"=-1), t0=31, T=45,
                       trace=TRUE, simplify=TRUE)
  , "Simulation has ended @t = 45"
)

## check construction and selection from "simEpidataCSlist"
mysim_from_list <- mysims[[1]]
invisible(capture.output(
    mysim_single <- eval(replace(attr(mysims, "call"), "nsim", 1))
))
mysim_from_list$runtime <- mysim_single$runtime <- NULL
expect_equivalent(mysim_single, mysim_from_list)

## check equivalence of Lambdag from simulation and residuals via twinstim
expect_equal(
    suppressMessages(surveillance:::residuals.twinstim(surveillance:::as.twinstim.simEpidataCS(mysims[[1]]))),
    residuals(mysims[[1]])
)
