context("Spatial interaction functions")


### test bundle

myexpectation <- function (siaf, intrfr, intrderivr, pargrid, type = 1, ...)
{
    ## check analytical intrfr specification against numerical approximation
    if (!missing(intrfr)) apply(pargrid, 1, function (pars) expect_warning(
        polyCub::checkintrfr(intrfr, siaf$f, pars, type, center=c(0,0),
                              rs=c(1,2,5,10,20,50)),
        NA, label = "polyCub::checkintrfr()"))

    ## also check intrfr for deriv
    if (!missing(intrderivr)) for (paridx in seq_along(intrderivr))
        apply(pargrid, 1, function (pars) expect_warning(
            polyCub::checkintrfr(intrderivr[[paridx]],
                                  function (...) siaf$deriv(...)[,paridx],
                                  pars, type, center=c(0,0),
                                  rs=c(1,2,5,10,20,50)),
            NA,
            label = paste0("polyCub::checkintrfr() for deriv[,",paridx,"]")))

    ## check deriv, F, Deriv against numerical approximations
    checksiafres <- surveillance:::checksiaf(siaf, pargrid, type, ...)
    for (i in which(!sapply(checksiafres, is.null)))
        expect_true(unique(attr(checksiafres[[i]], "all.equal")),
                    label=names(checksiafres)[i])
}


### test all pre-defined spatial interaction functions

test_that("Gaussian implementation agrees with numerical approximation",
          myexpectation(siaf.gaussian(F.adaptive=TRUE),
                        pargrid=t(log(0.5)),
                        tolerance=0.0005, method="midpoint", dimyx=250))

test_that("Power-law implementation agrees with numerical approximation",
          myexpectation(siaf.powerlaw(engine = "R"),
                        surveillance:::intrfr.powerlaw,
                        list(surveillance:::intrfr.powerlaw.dlogsigma,
                             surveillance:::intrfr.powerlaw.dlogd),
                        pargrid=cbind(0.5,log(c(0.1,1,2))),
                        tolerance=0.0005, method="SV", nGQ=13))

test_that("Lagged power-law implementation agrees with numeric results",
          myexpectation(siaf.powerlawL(engine = "R"),
                        surveillance:::intrfr.powerlawL,
                        list(surveillance:::intrfr.powerlawL.dlogsigma,
                             surveillance:::intrfr.powerlawL.dlogd),
                        pargrid=cbind(-0.5,log(c(0.1,1,2))),
                        tolerance=0.0005, method="midpoint", dimyx=250))

test_that("Student implementation agrees with numerical approximation",
          myexpectation(siaf.student(engine = "R"),
                        surveillance:::intrfr.student,
                        list(surveillance:::intrfr.student.dlogsigma,
                             surveillance:::intrfr.student.dlogd),
                        pargrid=cbind(0.5,log(c(0.1,1,2))),
                        tolerance=0.0005, method="SV", nGQ=5))

test_that("Step kernel implementation agrees with numerical approximation",
          myexpectation(siaf.step(c(0.1,0.5,1)),
                        pargrid=-t(c(0.5,0.1,0.2)),
                        tolerance=0.0005, method="midpoint", dimyx=350))


## ## plot the polygon on which F and Deriv are tested (to choose parameters)
## showsiaf <- function (siaf, pars) {
##     data("letterR", package="spatstat", envir=environment())
##     poly <- spatstat::shift.owin(letterR, -c(3,2))
##     plotpolyf(poly, siaf$f, pars, print.args=list(split=c(1,1,2,1), more=TRUE))
##     plotpolyf(poly, function (...) siaf$deriv(...)[,1], pars, print.args=list(split=c(2,1,2,1)))
## }
## showsiaf(siaf.student(), c(0.5,-0.5))


### test new C-implementations of F and Deriv functions

expect_equal_CnR <- function (siafgen, pargrid)
{
    polydomain <- spatstat::shift.owin(spatstat::letterR, -c(3,2))
    siafR <- siafgen(engine = "R")
    siafC <- siafgen(engine = "C")
    ## check F
    resF <- apply(pargrid, 1, function (pars)
        c(C = siafC$F(polydomain, , pars),
          R = siafR$F(polydomain, , pars)))
    expect_equal(object = resF["C",], expected = resF["R",],
                 label = "C-version of F", expected.label = "R-version of F")
    ## check Deriv
    resDeriv <- apply(pargrid, 1, function (pars)
        c(siafC$Deriv(polydomain, , pars),
          siafR$Deriv(polydomain, , pars)))
    p <- siafR$npars
    expect_equal(object = resDeriv[seq_len(p),], expected = resDeriv[p+seq_len(p),],
                 label = "C-version of Deriv", expected.label = "R-version of Deriv")
}

test_that("siaf.powerlaw() engines agree", {
    expect_equal_CnR(siafgen = siaf.powerlaw,
                     pargrid = cbind(0.5,log(c(0.1,1,2))))
})

test_that("siaf.student() engines agree", {
    expect_equal_CnR(siafgen = siaf.student,
                     pargrid = cbind(0.5,log(c(0.1,1,2))))
})

test_that("siaf.powerlawL() engines agree", {
    expect_equal_CnR(siafgen = siaf.powerlawL,
                     pargrid = cbind(-0.5,log(c(0.1,1,2))))
})
