context("earsC method")

test_that("earsC returns a sts object", {
  #Sim data and convert to sts object
  disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208, A = 1,
                                alpha = 1, beta = 0, phi = 0,
                                frequency = 1, state = NULL, K = 1.7)
  stsObj = disProg2sts( disProgObj)
  
  res1 <- earsC(stsObj, control = list(range = 20:208, method = "C1"))
  res2 <- earsC(stsObj, control = list(range = 20:208, method = "C2",
                                       alpha = 0.05))

  res3 <- earsC(stsObj, control = list(range = 20:208, method = "C3", sigma = 0.5))
  
  expect_is(res1, "sts")
  expect_is(res2, "sts")
  expect_is(res3, "sts")
  
  data("salmNewport")
  in2011 <- which(isoWeekYear(epoch(salmNewport))$ISOYear == 2011)
  salmNewportGermany <- aggregate(salmNewport, by = "unit")
  control <- list(range = in2011, method = "C1", alpha = 0.05)
  surv <- earsC(salmNewportGermany, control = control)
  
  expect_is(surv, "sts")
  expect_true(max(surv@upperbound[1:4] -
                  c(3.278854, 3.278854, 3.436517, 3.855617)) < 0.000001)
})

test_that("earsC returns error messages",{
  data("salmNewport")
  salmNewportGermany <- aggregate(salmNewport, by = "unit")
  control <- list(range = length(salmNewportGermany), method = "C1", alpha = 0.05,
                  baseline = 2)
  expect_error(earsC(salmNewportGermany, control = control), 
               "Minimum baseline to use is 3.")
  
  control <- list(range = length(salmNewportGermany), method = "C1", alpha = 0.05,
                  minSigma = - 2)
  expect_error(earsC(salmNewportGermany, control = control), 
               "The minimum sigma parameter")
  
  in2011 <- which(isoWeekYear(epoch(salmNewport))$ISOYear == 2011)
  control <- list(range = in2011, method = "C1", alpha = 0.05,
                  baseline = 1500)
  expect_error(earsC(salmNewportGermany, control = control), 
               "The vector of observed is too short!")
  
})

test_that("The range is well defined",{
  data("salmNewport")
  salmNewportGermany <- aggregate(salmNewport, by = "unit")
  control <- list(range = length(salmNewportGermany), method = "C1", alpha = 0.05,
                  baseline = 2)
  surv <- earsC(salmNewportGermany, 
                control = list(method = "C1",
                               baseline = 10))
  
  expect_true(length(surv@upperbound) == length(salmNewportGermany@observed) - 10)
  
})
