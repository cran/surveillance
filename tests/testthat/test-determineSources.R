data("imdepi")

test_that("determineSources() yields same result as old implementation", {
    sources0 <- surveillance:::determineSources.epidataCS(imdepi, method = "R")
    expect_identical(sources0, imdepi$events$.sources)
    sources1 <- surveillance:::determineSources(
        imdepi$events$time, imdepi$events$eps.t,
        coordinates(imdepi$events), imdepi$events$eps.s,
        imdepi$events$type, imdepi$qmatrix
        )
    expect_identical(sources1, imdepi$events$.sources)
    sources2 <- surveillance:::determineSources.epidataCS(imdepi, method = "C")
    expect_identical(sources2, imdepi$events$.sources)
})
