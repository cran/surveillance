context("Determine list of potential sources in \"epidataCS\"")

data("imdepi")

test_that("determineSourcesC() yields same result as old implementation", {
    sources0 <- determineSources.epidataCS(imdepi, method = "R")
    expect_identical(sources0, imdepi$events$.sources)
    sources1 <- determineSources(imdepi$events$time, imdepi$events$eps.t,
                                 coordinates(imdepi$events), imdepi$events$eps.s,
                                 imdepi$events$type, imdepi$qmatrix)
    expect_identical(sources1, imdepi$events$.sources)
    sources2 <- determineSources.epidataCS(imdepi, method = "C")
    expect_identical(sources2, imdepi$events$.sources)
})
