data("salmonella.agona")
# sts object
lala <- paste(salmonella.agona$start[1],salmonella.agona$start[2],"1",sep=" ")
firstMonday <- as.POSIXlt(lala, format = "%Y %W %u")
salm.ts <- salmonella.agona$observed
dates <- as.Date(firstMonday) + 7 * 0:(length(salm.ts) - 1)
start=c(salmonella.agona$start[1],salmonella.agona$start[2])
salm <- new("sts",epoch = as.numeric(dates), start = start, freq = 52,
observed = salm.ts, epochAsDate = TRUE)


###
## WEIGHTS FUNCTION
###

test_that("gamma = 1 if everything below the threshold",{
  s <- rep(0,10)
  weightsThreshold <- 0
  weights <- algo.farrington.assign.weights(s,weightsThreshold)
  expect_equal(weights,rep(1,10))
})

test_that(" A case that was checked by hand",{
  s <- rep(2,10)
  s[1:5] <- 0
  weightsThreshold <- 0
  weights <- algo.farrington.assign.weights(s,weightsThreshold)
  expect_equal(weights[1:5],rep(1.6,5))
  expect_equal(weights[6:10],rep(0.4,5))
})


###
## RESIDUALS FUNCTION
###

test_that(" residuals should be zero",{
  x <- rpois(10,1)
  y <- exp(x)
  model <- glm(y~x,family = quasipoisson(link="log"))
  phi <- max(summary(model)$dispersion,1)
  s <- anscombe.residuals(model,phi)
  expect_equal(as.numeric(s),rep(0,10))
})

test_that(" residuals should not be zero",{
  x <- rpois(1000,1)
  y <- exp(x)+runif(1)
  model <- glm(y~x,family = quasipoisson(link="log"))
  phi <- max(summary(model)$dispersion,1)
  s <- anscombe.residuals(model,phi)
  expect_true(mean(s)>0)
})


###
## FORMULA FUNCTION
###

test_that("We get the right formula",{
  expect_identical(surveillance:::formulaGLM(populationOffset=FALSE,timeBool=TRUE,factorsBool=FALSE),
                   "response ~ 1+wtime")
  expect_identical(surveillance:::formulaGLM(populationOffset=FALSE,timeBool=FALSE,factorsBool=FALSE),
                   "response ~ 1")
  expect_identical(surveillance:::formulaGLM(populationOffset=TRUE,timeBool=TRUE,factorsBool=FALSE),
                   "response ~ 1+wtime+offset(log(population))")
  expect_identical(surveillance:::formulaGLM(populationOffset=TRUE,timeBool=TRUE,factorsBool=TRUE),
                   "response ~ 1+wtime+offset(log(population))+seasgroups")
})


###
## REFERENCE TIME POINTS FUNCTION
###

test_that("We get the expected timepoints with weekly data",{
  # Case with weekly data with dates
  dayToConsider <- as.Date("2013-06-06")
  b <- 3
  freq <- 52
  epochAsDate <- TRUE
  epochStr <- "week"
  lala <- surveillance:::algo.farrington.referencetimepoints(dayToConsider,b=b,freq=freq,epochAsDate,epochStr)
  # Do we get the same day as dayToConsider?
  expect_equal(as.numeric(format(lala, "%w")),rep(4,4))
  # Actually for this example I know the dates one should get
  expect_equal(sort(lala),sort(c(as.Date("2010-06-03"),as.Date("2013-06-06"),as.Date("2012-06-07"),as.Date("2011-06-09"))))
})

test_that("We get the expected timepoints with monthly data",{
  dayToConsider <- 48
  b <- 3
  freq <- 12
  epochAsDate <- FALSE
  epochStr <- "month"
  lala <- surveillance:::algo.farrington.referencetimepoints(dayToConsider,b=b,freq=freq,epochAsDate,epochStr)
  expect_equal(lala,c(48,36,24,12))
})
test_that("We get an error when going too many years back",{
  dayToConsider <- 48
  b <- 3
  freq <- 12
  epochAsDate <- FALSE
  epochStr <- "month"
  expect_true(any(surveillance:::algo.farrington.referencetimepoints(dayToConsider,b=8,freq=freq,epochAsDate,epochStr) < 1))

  # apply code
   control1 <-  list(range=250,noPeriods=10,populationOffset=FALSE,
                     fitFun="algo.farrington.fitGLM.flexible",
                     b=10,w=3,weightsThreshold=2.58,
                     pastWeeksNotIncluded=26,
                     pThresholdTrend=1,trend=TRUE,
                     thresholdMethod="muan",alpha=0.05,glmWarnings=FALSE)
  expect_error(farringtonFlexible(salm,control=control1),"Some reference")
})


###
## FIT GLM FUNCTION
###

# Case with convergence
control<-  list(range=250,noPeriods=10,populationOffset=TRUE,
                fitFun="algo.farrington.fitGLM.flexible",
                b=40,w=3,weightsThreshold=2.58,
                pastWeeksNotIncluded=26,
                pThresholdTrend=1,trend=TRUE,
                thresholdMethod="muan",alpha=0.05,glmWarnings=FALSE)
 response=salm@observed[1:120]
dataGLM <- data.frame(response=response,wtime=1:120,
	                  population=runif(120)*100,
                      seasgroups=as.factor(rep(1:12,10)))

arguments <- list(dataGLM=dataGLM,
                   timeTrend=TRUE,
                   populationOffset=TRUE,
                   factorsBool=TRUE,reweight=TRUE,
                   weightsThreshold=0.5,glmWarnings=control$glmWarnings,
				   control=control)
model <- do.call(surveillance:::algo.farrington.fitGLM.flexible, args=arguments)

test_that("The fit glm function gives the right class of output",{
  expect_inherits(model, "glm")
})

test_that("The fit glm function gives as many coefficients as expected",{
  expect_equal(dim(summary(model)$coefficients)[1],
               length(levels(dataGLM$seasgroups))-1+1+1)
})

test_that("wtime, response, phi and weights were added to the model",{
  expect_false(is.null(model$phi))
  expect_false(is.null(model$wtime))
  expect_false(is.null(model$response))
  expect_false(is.null(model$population))
  expect_false(is.null(model$weights))
})

test_that("reweighting was done",{
  expect_true(all(model$weights!=1))
})

test_that("there are no weights if very high threshold",{
  arguments$reweight <- TRUE
  arguments$weightsThreshold <- 100000
  model <- do.call(surveillance:::algo.farrington.fitGLM.flexible, args=arguments)
  expect_true(all(model$weights==1))
})

test_that("there is not a too small overdispersion",{
  expect_true(model$phi>=1)
})


###
## BLOCKS FUNCTION
###

referenceTimePoints <- c(as.Date("2010-06-03"),as.Date("2013-06-06"),as.Date("2012-06-07"),as.Date("2011-06-09"))
firstDay <- as.Date("1990-06-07")
vectorOfDates <- dates <- as.Date(firstDay) + 7 * 0:1300
freq <- 52
dayToConsider <- as.Date("2013-06-06")
b <- 3
w <- 3
epochAsDate <- TRUE

# p=1
p <- 1
lala <- surveillance:::blocks(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,epochAsDate)
test_that("the reference window has the right length",{
  expect_equal(length(vectorOfDates[is.na(lala)==FALSE&lala==p]),w+1+b*(2*w+1))

  # p>1
  p <- 8
  lala <- surveillance:::blocks(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,epochAsDate)
  # reference windows
  expect_equal(length(vectorOfDates[is.na(lala)==FALSE&lala==p]),w+1+b*(2*w+1))
})

lili <- as.factor(lala[is.na(lala)==FALSE])

test_that("there are as many levels as expected",{
  expect_equal(length(levels(lili)),p)

})
p <- 8
lala <- surveillance:::blocks(referenceTimePoints,vectorOfDates,freq,dayToConsider,b,w,p,epochAsDate)
lili <- as.factor(lala[is.na(lala)==FALSE])
lolo <- lili[lili!=p]
test_that("periods of roughly the same length each year",{
    expect_equal(as.numeric(abs(diff(table(lolo))[1:(p-2)])<=b),rep(1,(p-2)))
})


###
## THRESHOLD FUNCTION FARRINGTON
###

predFit <- 5
predSeFit <- 0.2
wtime <- 380
skewness.transform <- "2/88"
alpha <- 0.05
y <- 8
method <- "delta"
phi <- 1

test_that("the function recognizes wrong exponents",{
  expect_error(surveillance:::algo.farrington.threshold.farrington(
    predFit, predSeFit, phi, skewness.transform, alpha, y, method
  ), "proper exponent")
})

test_that("some results we know are found",{
  skewness.transform <- "none"
  lala <- surveillance:::algo.farrington.threshold.farrington(
      predFit, predSeFit, phi, skewness.transform, alpha, y, method
  )
  # Should always be ok
  lala <- as.numeric(lala)
  expect_true(lala[3]<=1&lala[1]>=0)
  expect_true(lala[2]>lala[1])
  expect_true(lala[1]>=0)

  # Here we know the results
  expect_equal(as.numeric(lala), c(1.3073128, 8.6926872, 0.0907246, 0.8124165),
               tolerance = 1e-6, scale = 1)

  # Here we calculated some examples
  skewness.transform <- "1/2"
  lala <- surveillance:::algo.farrington.threshold.farrington(
    predFit, predSeFit, phi, skewness.transform, alpha, y, method
  )
  expect_equal(as.numeric(lala), c(1.9891097, 9.3744842, 0.1189986, 0.6857951),
               tolerance = 1e-6, scale = 1)

  skewness.transform <- "2/3"
  lala <- surveillance:::algo.farrington.threshold.farrington(
    predFit, predSeFit, phi, skewness.transform, alpha, y, method
  )
  expect_equal(as.numeric(lala), c(1.8084477, 9.1154825, 0.1094727, 0.7289546),
               tolerance = 1e-6, scale = 1)
})


###
## THRESHOLD FUNCTION NOUFAILY
###

predFit <- log(5)
predSeFit <- log(2)
wtime <- 380
skewness.transform <- "none"
alpha <- 0.05
y <- 11
phi <- 1.5
method <- "muan"
lala <- surveillance:::algo.farrington.threshold.noufaily(
  predFit, predSeFit, phi, skewness.transform, alpha, y, method
)
test_that("some results we know are found",{
  # Should always be ok
  lala <- as.numeric(lala)
  expect_true(lala[3]<=1&lala[1]>=0)
  expect_true(lala[2]>lala[1])
  expect_true(lala[1]>=0)

  # Here we calculated some examples
  expect_equal(as.numeric(lala), c(8.0000000, 24.0000000, 0.8597797, 0.4193982),
               tolerance = 1e-6, scale = 1)

  phi <- 1.0
  method <- "muan"
  lala <- surveillance:::algo.farrington.threshold.noufaily(
    predFit, predSeFit, phi, skewness.transform, alpha, y, method
  )
  expect_equal(as.numeric(lala), c(9.0000000, 22.0000000, 0.9093099, 0.4605347),
               tolerance = 1e-6, scale = 1)

  phi <- 1.5
  method <- "nbPlugin"
  lala <- surveillance:::algo.farrington.threshold.noufaily(
    predFit, predSeFit, phi, skewness.transform, alpha, y, method
  )
  expect_equal(as.numeric(lala), c(1.00000000, 10.00000000, 0.03763657,  1.11918153),
               tolerance = 1e-6, scale = 1)

  phi <- 1.0
  method <- "nbPlugin"
  lala <- surveillance:::algo.farrington.threshold.noufaily(
    predFit, predSeFit, phi, skewness.transform, alpha, y, method
  )
  expect_equal(as.numeric(lala), c(2.00000000, 9.00000000, 0.01369527, 1.27061541),
               tolerance = 1e-6, scale = 1)
})


###
## DATA GLM FUNCTION
###

b <- 3
freq <- 52
dayToConsider <- as.Date("2013-05-30")
epochAsDate <- TRUE
epochStr <- "week"
firstDay <- as.Date("1990-06-07")
vectorOfDates <- dates <- as.Date(firstDay) + 7 * 0:1300
w <- 3
noPeriods <- 10
observed <- rnorm(1301)+runif(1301)+30
population <- rnorm(1301)+10
verbose <- FALSE
pastWeeksNotIncluded <- w
k <- 1200

lala <- surveillance:::algo.farrington.data.glm(dayToConsider, b, freq,
                                 epochAsDate,epochStr,
                                 vectorOfDates,w,noPeriods,
                                 observed,population,
                                 verbose,pastWeeksNotIncluded,k)

test_that("the output is a data.frame",{
  expect_inherits(lala, "data.frame")
})

test_that("the data frame contains all variables",{
  expect_identical(names(lala), c("response", "wtime","population","seasgroups","vectorOfDates"))
})

test_that("the time variable is ok with diff 1",{
  expect_equal(diff(lala$wtime), rep(1,length(lala$wtime)-1))
})

test_that("the factor variable has the right number of levels",{
  expect_equal(nlevels(lala$seasgroups), noPeriods)
})

observed[1150] <- NA
lala <- surveillance:::algo.farrington.data.glm(dayToConsider, b, freq,
                                 epochAsDate,epochStr,
                                 vectorOfDates,w,noPeriods,
                                 observed,population,
                                 verbose,pastWeeksNotIncluded,k)

test_that("the data frame has the right dimensions",{
  expect_equal(dim(lala),c(156,5))
})


###
## GLM FUNCTION
###

dataGLM <- lala
timeTrend <- TRUE
populationOffset <- TRUE
factorsBool <- TRUE
reweight <- TRUE
weightsThreshold <- 1
pThresholdTrend <- 1
b <- 3
noPeriods <- 10
typePred <- "link"
fitFun <- "algo.farrington.fitGLM.flexible"
glmWarnings <- FALSE
epochAsDate <- TRUE
dayToConsider <- as.Date("2013-05-30")
diffDates <- 7
populationNow <- 10

test_that("the output has the needed variables",{
  finalModel <- surveillance:::algo.farrington.glm(dataGLM,timeTrend,populationOffset,factorsBool,
                                    reweight,weightsThreshold,pThresholdTrend,b,
                                    noPeriods,typePred,fitFun,glmWarnings,epochAsDate,
                                    dayToConsider,diffDates,populationNow,verbose=FALSE)
  expect_identical(names(finalModel), c("pred","doTrend","coeffTime","phi"))
})

test_that("no time trend in no time trend",{
  pThresholdTrend <- 1
  b <- 2
  finalModel <- surveillance:::algo.farrington.glm(dataGLM,timeTrend,populationOffset,factorsBool,
                                    reweight,weightsThreshold,pThresholdTrend,b,
                                    noPeriods,typePred,fitFun,glmWarnings,epochAsDate,
                                    dayToConsider,diffDates,populationNow,verbose=FALSE)
  expect_false(finalModel$doTrend)
})


###
## ALARMS
###

test <- farringtonFlexible(salm,control=list(thresholdMethod="nbPlugin",alpha=0.1))
test_that("there are only alarms when expected",{
  # No alarm when observed is 0
  expect_true(sum(test@alarm[test@observed==0])==0)
  # No alarm when the observed counts are UNDER the threshold
  expect_true(sum(observed(test)>upperbound(test),na.rm=TRUE)==sum(test@alarm==TRUE))
})


###
## NO CONVERGENCE
###

timeSeries <- rep(0,698)
timeSeries[696] <- 1

algoControl <- list(noPeriods=10,alpha = 0.01,verbose = F,
                    b=5,w=4,weightsThreshold=2.58,pastWeeksNotIncluded=26,
                    pThresholdTrend=1,thresholdMethod='nbPlugin',limit54 = c(4,5),
                    range = (length(timeSeries) - 1):length(timeSeries), glmWarnings = FALSE)
seriesSTSObject <- new('sts', observed = timeSeries,
                       epoch = as.numeric(seq(as.Date('2001-01-01'),length.out=length(timeSeries), by='1 week')),
                       epochAsDate = TRUE)
test_that("The code does not produce any error",{
# It is ok if the code does not produce any error
 expect_warning(farringtonFlexible(seriesSTSObject, control = algoControl))
})


###
## NA
###

timeSeries <- rnorm(698)*10+runif(698)*100+30

algoControl <- list(noPeriods=10,alpha = 0.01,verbose = F,
                    b=5,w=4,weightsThreshold=2.58,pastWeeksNotIncluded=NULL,
                    pThresholdTrend=1,thresholdMethod='nbPlugin',limit54 = c(4,5),
                    range = (length(timeSeries) - 1):length(timeSeries), glmWarnings = FALSE)
seriesSTSObject <- sts(timeSeries,
                       epoch = seq(as.Date('2001-01-01'), length.out=length(timeSeries), by='1 week'))
test_that("The code does not produce any error",{
  farringtonFlexible(seriesSTSObject, control = algoControl)

  results1 <- farringtonFlexible(seriesSTSObject, control = algoControl)
  expect_inherits(results1, "sts")
  seriesSTSObject@observed[680:690] <- NA
  results2 <- farringtonFlexible(seriesSTSObject, control = algoControl)
  expect_inherits(results2, "sts")
})


###
## multivariate example with populationOffset
###

msts <- sts(cbind(series1 = timeSeries, series2 = timeSeries),
            population = c(1,2)*10000)
algoControl$populationOffset <- TRUE
out <- farringtonFlexible(msts, control = algoControl)
test_that("population offsets of multivariate sts are handled correctly", {
  ## resulting upperbounds must be the same as different offsets just rescale the intercept
  expect_equal(out@upperbound[,1], out@upperbound[,2])
})
## failed in surveillance <= 1.19.1 (used series 1 offset for both fits)
