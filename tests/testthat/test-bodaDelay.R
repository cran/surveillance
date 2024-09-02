###
## Checking the provided reporting triangle
###

data('salmAllOnset')

# Control slot for the proposed algorithm with D=10 correction
rangeTest <- 410:412
alpha <- 0.05
controlDelay <-  list(range = rangeTest, b = 4, w = 3,
                      pastAberrations = TRUE, mc.munu=10, mc.y=10,
                      verbose = FALSE,populationOffset=FALSE,
                      alpha = alpha, trend = TRUE,
                      limit54=c(0,50),
                      noPeriods = 10, pastWeeksNotIncluded = 26,
                      delay=TRUE)
test_that("The absence of reporting triangle throws an error",{
  data("salmNewport")
  expect_error(bodaDelay(salmNewport, controlDelay),"You have to")
})
test_that("The function spots incorrect reporting triangles",{
  stsFake <- salmAllOnset
  stsFake@control$reportingTriangle$n <- head(stsFake@control$reportingTriangle$n,n=10)
  expect_error(bodaDelay(stsFake, controlDelay),"The reporting triangle number")
  stsFake <- salmAllOnset
  stsFake@control$reportingTriangle$n[1,] <- stsFake@control$reportingTriangle$n[1,]/2
  expect_error(bodaDelay(stsFake, controlDelay),"The reporting triangle is wrong")
})


###
## Data glm function
###

epochAsDate <- TRUE
epochStr <- "week"
freq <- 52
b <- controlDelay$b
w <- controlDelay$w
populationOffset <- controlDelay$populationOffset
noPeriods <- controlDelay$noPeriods
verbose <- controlDelay$verbose
reportingTriangle <- salmAllOnset@control$reportingTriangle
timeTrend <- controlDelay$trend
alpha <- controlDelay$alpha
populationOffset <- controlDelay$populationOffset
factorsBool <- controlDelay$factorsBool
pastAberrations <- controlDelay$pastAberrations
glmWarnings <- controlDelay$glmWarnings
delay <- controlDelay$delay
k <- controlDelay$k
verbose <- controlDelay$verbose
pastWeeksNotIncluded <- controlDelay$pastWeeksNotIncluded
mc.munu <- controlDelay$mc.munu
mc.y <- controlDelay$mc.y
vectorOfDates <- as.Date(salmAllOnset@epoch, origin="1970-01-01")
dayToConsider <- vectorOfDates[rangeTest[1]]
observed <- salmAllOnset@observed
population <- salmAllOnset@populationFrac

dataGLM <- surveillance:::bodaDelay.data.glm(dayToConsider=dayToConsider,
                              b=b, freq=freq,
                              epochAsDate=epochAsDate,
                              epochStr=epochStr,
                              vectorOfDates=vectorOfDates,w=w,
                              noPeriods=noPeriods,
                              observed=observed,population=population,
                              verbose=verbose,
                              pastWeeksNotIncluded=pastWeeksNotIncluded,
                              reportingTriangle=reportingTriangle,
                              delay=delay)
delay <- FALSE
dataGLMNoDelay <- surveillance:::bodaDelay.data.glm(dayToConsider=dayToConsider,
                                  b=b, freq=freq,
                                  epochAsDate=epochAsDate,
                                  epochStr=epochStr,
                                  vectorOfDates=vectorOfDates,w=w,
                                  noPeriods=noPeriods,
                                  observed=observed,population=population,
                                  verbose=verbose,
                                  pastWeeksNotIncluded=pastWeeksNotIncluded,
                                  reportingTriangle=reportingTriangle,
                                  delay=delay)

test_that("the output is a data.frame",{
  expect_inherits(dataGLM, "data.frame")
  expect_inherits(dataGLMNoDelay, "data.frame")
})

test_that("the data frame contains all variables",{
  expect_identical(names(dataGLM), c("response", "wtime","population","seasgroups","vectorOfDates","delay"))
  expect_identical(names(dataGLMNoDelay), c("response", "wtime","population","seasgroups","vectorOfDates"))
  })

test_that("the variables have the right class",{
  expect_inherits(dataGLM$response, "numeric")
  expect_inherits(dataGLM$wtime, "numeric")
  expect_inherits(dataGLM$population, "numeric")
  expect_inherits(dataGLM$seasgroups, "factor")
  expect_inherits(dataGLM$vectorOfDates, "Date")
  expect_inherits(dataGLM$delay, "numeric")

  expect_inherits(dataGLMNoDelay$response, "numeric")
  expect_inherits(dataGLMNoDelay$wtime, "numeric")
  expect_inherits(dataGLMNoDelay$population, "numeric")
  expect_inherits(dataGLMNoDelay$seasgroups, "factor")
  expect_inherits(dataGLMNoDelay$vectorOfDates, "Date")
})

test_that("the time variable is ok with diff 1",{
  delayWtime <- as.numeric(levels(as.factor(dataGLM$wtime)))
  expect_equal(diff(delayWtime), rep(1,length(delayWtime)-1))
  expect_equal(diff(dataGLMNoDelay$wtime), rep(1,length(dataGLMNoDelay$wtime)-1))
})

test_that("the factor variable has the right number of levels",{
  expect_equal(nlevels(dataGLM$seasgroups), noPeriods)
  expect_equal(nlevels(dataGLMNoDelay$seasgroups), noPeriods)
})


###
## Fit glm function
###

argumentsGLM <- list(dataGLM=dataGLM,reportingTriangle=reportingTriangle,
                     timeTrend=timeTrend,alpha=alpha,
                     populationOffset=populationOffset,
                     factorsBool=TRUE,pastAberrations=FALSE,
                     glmWarnings=glmWarnings,
                     verbose=verbose,delay=delay,k=k,control=controlDelay)

if(surveillance.options("allExamples") && require("INLA")) { # needs to be attached
  argumentsGLM$inferenceMethod <- "INLA"
  model <- do.call(surveillance:::bodaDelay.fitGLM, args=argumentsGLM)
  test_that("the fitGLM function gives the right class of output",{
    expect_inherits(model, "inla")
  })
}

argumentsGLM$inferenceMethod <- "asym"
model <- do.call(surveillance:::bodaDelay.fitGLM, args=argumentsGLM)
test_that("the fitGLM function gives the right class of output",{
  expect_inherits(model, "negbin")
})


###
## formula function
###

test_that("We get the right formula",{
  expect_identical(surveillance:::formulaGLMDelay(timeBool=TRUE,factorsBool=FALSE),
                   "response ~ 1+wtime")
  expect_identical(surveillance:::formulaGLMDelay(timeBool=FALSE,factorsBool=FALSE),
                   "response ~ 1")
  expect_identical(surveillance:::formulaGLMDelay(timeBool=TRUE,factorsBool=FALSE),
                   "response ~ 1+wtime")
  expect_identical(surveillance:::formulaGLMDelay(timeBool=TRUE,factorsBool=TRUE),
                   "response ~ 1+wtime+as.factor(seasgroups)")
  expect_identical(surveillance:::formulaGLMDelay(timeBool=TRUE,factorsBool=TRUE,delay=TRUE),
                   "response ~ 1+wtime+as.factor(seasgroups)+as.factor(delay)")
  expect_identical(surveillance:::formulaGLMDelay(timeBool=TRUE,factorsBool=FALSE,outbreak=TRUE),
                   "response ~ 1+wtime+f(outbreakOrNot,model='linear', prec.linear = 1)")
})
