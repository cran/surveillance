## ----SETUP, include = FALSE--------------------------------------------------------
## create directories for plots and cache
dir.create("plots", showWarnings=FALSE)
dir.create("monitoringCounts-cache", showWarnings=FALSE)
## load packages
library('surveillance')
library('gamlss')

## ----echo=FALSE--------------------------------------------------------------------
data("salmNewport")

## ----echo=FALSE--------------------------------------------------------------------
stopifnot(
all.equal(observed(salmNewport),
          observed(as(as(salmNewport, "ts"), "sts")))
)

## ----echo=FALSE--------------------------------------------------------------------
# This code is the one used for the Salmon et al. (2016) JSS article.
# Using this code all examples from the article can be reproduced.
# computeALL is FALSE to avoid the computationally intensive parts
# of the code (simulations to find a threshold value for categoricalCUSUM,
# INLA-driven BODA) but one can set it to TRUE to have it run.
computeALL <- FALSE

## ----NewportPlot-simple, fig.keep = 'none'-----------------------------------------
plot(salmNewport, type = observed ~ time,
     xaxis.tickFreq = list("%m" = atChange, "%G" = atChange),
     xaxis.labelFreq = list("%Y" = atMedian), xaxis.labelFormat = "%Y")

## ----unitPlot-simple, echo = FALSE, fig.keep = 'none'------------------------------
plot(salmNewport, units = 2:3)

## ----EARS, fig.keep='none'---------------------------------------------------------
in2011 <- which(isoWeekYear(epoch(salmNewport))$ISOYear == 2011)
salmNewportGermany <- aggregate(salmNewport, by = "unit")
control <- list(range = in2011, method = "C1", alpha = 0.05)
surv <- earsC(salmNewportGermany, control = control)
plot(surv)

## ----farHead-----------------------------------------------------------------------
con.farrington <- list(
    range = in2011, noPeriods = 1,
    b = 4, w = 3, weightsThreshold = 1,
    pastWeeksNotIncluded = 3, pThresholdTrend = 0.05,
    thresholdMethod = "delta"
)
con.noufaily   <- list(
    range = in2011, noPeriods = 10,
    b = 4, w = 3, weightsThreshold = 2.58,
    pastWeeksNotIncluded = 26, pThresholdTrend = 1,
    thresholdMethod = "nbPlugin"
)

## ----echo=F------------------------------------------------------------------------
con.farrington$limit54 <- con.noufaily$limit54 <- c(0,50)  # for the figure

## ----------------------------------------------------------------------------------
salm.farrington <- farringtonFlexible(salmNewportGermany, con.farrington)
salm.noufaily   <- farringtonFlexible(salmNewportGermany, con.noufaily)

## ----farPlot-simple, echo = FALSE, fig.keep = 'none'-------------------------------
par(mfrow = c(1,2))
plot(salm.farrington)
plot(salm.noufaily)

## ----campyDE-simple, fig.keep='none'-----------------------------------------------
data("campyDE")
cam.sts <- sts(epoch = campyDE$date,
               observed = campyDE$case, state = campyDE$state)

plot(cam.sts, col = "mediumblue")
lines(campyDE$hum * 50, col = "white", lwd = 2)
axis(4, at = seq(0, 2500, by = 500), labels = seq(0, 50, by = 10))

## ----boda-cache, echo = FALSE, results='hide'--------------------------------------
if (computeALL) {
## The original results were produced using version 0.0-1458166556,
## and version 0.0-1485844051 from 2017-01-31 also worked. However:
## hoehle 2018-07-18: changed to prior="iid" as "rw1" crashes INLA >= 17.06.20.
## smeyer 2025-06-24: restored prior="rw1", working again with INLA 25.06.07.
library("INLA")
rangeBoda <- which(epoch(cam.sts) >= as.Date("2007-01-01"))
control.boda <- list(range = rangeBoda, X = NULL, trend = TRUE,
                     season = TRUE, prior = "rw1", alpha = 0.025,
                     mc.munu = 10000, mc.y = 1000,
                     samplingMethod = "marginals")
boda <- boda(cam.sts, control = control.boda)
save(list = c("boda", "control.boda", "rangeBoda"),
     file = "monitoringCounts-cache/boda.RData")
} else {
  load("monitoringCounts-cache/boda.RData")
}

## ----boda2-cache, echo = FALSE, results='hide'-------------------------------------
if (computeALL) {
covarNames <- c("l1.hum", "l2.hum", "l3.hum", "l4.hum",
                "newyears", "christmas", "O104period")
control.boda2 <- modifyList(control.boda,
                            list(X = campyDE[, covarNames], season = FALSE))
boda.covars <- boda(cam.sts, control = control.boda2)
save(list = c("boda.covars", "covarNames", "control.boda2"),
     file = "monitoringCounts-cache/boda.covars.RData")
} else {
  load("monitoringCounts-cache/boda.covars.RData")
}

## ----bPlot-simple, echo = FALSE, fig.keep = 'none'---------------------------------
plot(boda.covars)

## ----boda3, echo = FALSE-----------------------------------------------------------
control.far <- list(range=rangeBoda,b=4,w=5,alpha=0.025*2)
far <- farrington(cam.sts,control=control.far)
#Both farringtonFlexible and algo.bayes uses a one-sided interval just as boda.
control.far2 <-modifyList(control.far,list(alpha=0.025))
farflex <- farringtonFlexible(cam.sts,control=control.far2)
bayes <- suppressWarnings(bayes(cam.sts,control=control.far2))

## ----boda4, echo = FALSE-----------------------------------------------------------
# Small helper function to combine several equally long univariate sts objects
combineSTS <- function(stsList) {
 epoch <- as.numeric(epoch(stsList[[1]]))
 observed <- NULL
 alarm <- NULL
 for (i in 1:length(stsList)) {
   observed <- cbind(observed,observed(stsList[[i]]))
   alarm <- cbind(alarm,alarms(stsList[[i]]))
 }
 colnames(observed) <- colnames(alarm) <- names(stsList)
 res <- sts(epoch=as.numeric(epoch), epochAsDate=TRUE,
            observed=observed, alarm=alarm)
 return(res)
}

## ----alarmplot, fig.width=8, fig.height=4, out.width="\\linewidth", echo=FALSE-----
# Make an artificial object containing two columns - one with the boda output
# and one with the farrington output
cam.surv <- combineSTS(list(boda.covars=boda.covars,boda=boda,bayes=bayes,
                            farrington=far,farringtonFlexible=farflex))
par(mar=c(4,8,2.1,2),family="Times")
plot(cam.surv,type = alarm ~ time,lvl=rep(1,ncol(cam.surv)),
     alarm.symbol=list(pch=17, col="red2", cex=1,lwd=3),
     cex.axis=1,xlab="Time (weeks)",cex.lab=1,xaxis.tickFreq=list("%m"=atChange,"%G"=atChange),xaxis.labelFreq=list("%G"=at2ndChange),
     xaxis.labelFormat="%G")

## ----glrnb, results='hide'---------------------------------------------------------
phase1 <- which(isoWeekYear(epoch(salmNewportGermany))$ISOYear < 2011)
phase2 <- in2011
control <- list(range = phase2, c.ARL = 4, theta = log(2), ret = "cases",
                mu0 = list(S = 1, trend = TRUE, refit = FALSE))
salmGlrnb <- glrnb(salmNewportGermany, control = control)

## ----cat---------------------------------------------------------------------------
data("salmHospitalized")
isoWeekYearData <- isoWeekYear(epoch(salmHospitalized))
dataBefore2013 <- which(isoWeekYearData$ISOYear < 2013)
data2013 <- which(isoWeekYearData$ISOYear == 2013)
dataEarly2014 <- which(isoWeekYearData$ISOYear == 2014
                       & isoWeekYearData$ISOWeek <= 4)
phase1 <- dataBefore2013
phase2 <- c(data2013, dataEarly2014)

salmHospitalized.df <- cbind(as.data.frame(salmHospitalized),
                             weekNumber = isoWeekYearData$ISOWeek)
names(salmHospitalized.df) <- c("y", "t", "state", "alarm", "upperbound", "n",
                                "freq", "epochInPeriod", "weekNumber")

## ----catbis, results='hide'--------------------------------------------------------
vars <- c( "y", "n", "t", "epochInPeriod", "weekNumber")
m.bbin <- gamlss(cbind(y, n-y) ~ 1 + t
                 + sin(2 * pi * epochInPeriod) + cos(2 * pi * epochInPeriod)
                 + sin(4 * pi * epochInPeriod) + cos(4 * pi * epochInPeriod)
                 + I(weekNumber == 1) + I(weekNumber == 2),
                 sigma.formula =~ 1,
                 family = BB(sigma.link = "log"),
                 data = salmHospitalized.df[phase1, vars])

## ----cat2--------------------------------------------------------------------------
R <- 2
h <- 2
pi0 <- predict(m.bbin, newdata = salmHospitalized.df[phase2, vars],
               type = "response")
pi1 <- plogis(qlogis(pi0) + log(R))
pi0m <- rbind(pi0, 1 - pi0)
pi1m <- rbind(pi1, 1 - pi1)

## ----cat2bis-----------------------------------------------------------------------
populationHosp <- unname(cbind(
  population(salmHospitalized),
  population(salmHospitalized)))
observedHosp <- cbind(
  "Yes" = as.vector(observed(salmHospitalized)),
  "No" = as.vector(population(salmHospitalized) - observed(salmHospitalized)))
salmHospitalized.multi <- sts(
  frequency = 52, start = c(2004, 1), epoch = epoch(salmHospitalized),
  observed = observedHosp, population = populationHosp,
  multinomialTS = TRUE)

## ----cat2ter-----------------------------------------------------------------------
dBB.cusum <- function(y, mu, sigma, size, log = FALSE) {
  dBB(if (is.matrix(y)) y[1,] else y,
      if (is.matrix(y)) mu[1,] else mu,
      sigma = sigma, bd = size, log = log)
}

## ----cat3--------------------------------------------------------------------------
controlCat <- list(range = phase2, h = 2, pi0 = pi0m, pi1 = pi1m,
                   ret = "cases", dfun = dBB.cusum)
salmHospitalizedCat <- categoricalCUSUM(
    salmHospitalized.multi, control = controlCat,
    sigma = exp(m.bbin$sigma.coefficients))

## ----------------------------------------------------------------------------------
h.grid <- seq(1, 10, by = 0.5)

## ----cath-cache, echo = FALSE, results='hide'--------------------------------------
if (computeALL) {
simone <- function(sts, h) {
  y <- rBB(length(phase2), mu = pi0m[1, , drop = FALSE],
           bd = population(sts)[phase2, ],
           sigma = exp(m.bbin$sigma.coefficients),
           fast = TRUE)
  observed(sts)[phase2, ] <- cbind(y, population(sts)[phase2, 1] - y)
  one.surv <- categoricalCUSUM(
      sts, control = modifyList(controlCat, list(h = h)),
      sigma = exp(m.bbin$sigma.coefficients))
  return(any(alarms(one.surv)[, 1]))
}
set.seed(123)
nSims <- 1000
pMC <- sapply(h.grid, function(h) {
  mean(replicate(nSims, simone(salmHospitalized.multi, h)))
})

pMarkovChain <- sapply(h.grid, function(h) {
  TA <- LRCUSUM.runlength(mu = pi0m[1,,drop = FALSE],
                          mu0 = pi0m[1,,drop = FALSE],
                          mu1 = pi1m[1,,drop = FALSE],
                          n = population(salmHospitalized.multi)[phase2, ],
                          h = h, dfun = dBB.cusum,
                          sigma = exp(m.bbin$sigma.coef))
  return(tail(TA$cdf, n = 1))
})
save(pMC, file = "monitoringCounts-cache/pMC.RData")
save(pMarkovChain, file = "monitoringCounts-cache/pMarkovChain.RData")
} else {
load("monitoringCounts-cache/pMC.RData")
load("monitoringCounts-cache/pMarkovChain.RData")
}

## ----catF-simple, echo = FALSE, fig.keep = 'none'----------------------------------
plot(salmHospitalizedCat[,1])

## ----catARL-simple, echo = FALSE, fig.keep = 'none'--------------------------------
matplot(h.grid, cbind(pMC, pMarkovChain), type="l", lty=1:2, col=1)
abline(h=0.1, lty=5)
legend("center", c("Monte Carlo","Markov chain"), lty=1:2, bty="n")

## ----rotaPlot-simple, fig.keep='none'----------------------------------------------
data("rotaBB")
plot(rotaBB)

## ----------------------------------------------------------------------------------
rotaBB.df <- as.data.frame(rotaBB)
X <- with(rotaBB.df, cbind(intercept = 1, epoch,
                           sin1 = sin(2 * pi * epochInPeriod),
                           cos1 = cos(2 * pi * epochInPeriod)))

phase1 <- epoch(rotaBB) < as.Date("2009-01-01")
phase2 <- !phase1

library("MGLM")
## MGLMreg automatically takes the last class as ref so we reorder
order <- c(2:5, 1); reorder <- c(5, 1:4)
m0 <- MGLMreg(as.matrix(rotaBB.df[phase1, order]) ~ -1 + X[phase1, ],
              dist = "MN")

## ----------------------------------------------------------------------------------
m1 <- m0
m1@coefficients[1, ] <- m0@coefficients[1, ] + log(2)

pi0 <- t(predict(m0, newdata = X[phase2, ])[, reorder])
pi1 <- t(predict(m1, newdata = X[phase2, ])[, reorder])

## ----CATCUSUM----------------------------------------------------------------------
dfun <- function(y, size, mu, log = FALSE) {
  dmultinom(x = y, size = size, prob = mu, log = log)
}

h <- 2  # threshold for the CUSUM statistic
control <- list(range = seq(nrow(rotaBB))[phase2], h = h, pi0 = pi0,
                pi1 = pi1, ret = "value", dfun = dfun)
surv <- categoricalCUSUM(rotaBB,control=control)

## ----include = FALSE---------------------------------------------------------------
alarmDates <- epoch(surv)[which(alarms(surv)[,1])]
format(alarmDates,"%b %Y")

## ----CATCUSUMMC,echo=FALSE,eval=FALSE----------------------------------------------
# #Number of MC samples
# nSamples <- 1e4
# 
# #Do MC
# simone.stop <- function(sts, control) {
#   phase2Times <- seq(nrow(sts))[phase2]
#   #Generate new phase2 data from the fitted in control model
#   y <- sapply(1:length(phase2Times), function(i) {
#     rmultinom(n=1, prob=pi0[,i],size=population(sts)[phase2Times[i],1])
#   })
#   observed(sts)[phase2Times,] <- t(y)
#   one.surv <- categoricalCUSUM(sts, control=control)
#   #compute P(S<=length(phase2))
#   return(any(alarms(one.surv)[,1]>0))
# }
# 
# set.seed(1233)
# rlMN <- replicate(nSamples, simone.stop(rotaBB, control=control))
# mean(rlMN)  # 0.5002

## ----------------------------------------------------------------------------------
m0.dm <- MGLMreg(as.matrix(rotaBB.df[phase1, 1:5]) ~ -1 + X[phase1, ],
                dist = "DM")
c(m0@AIC, m0.dm@AIC)

## ----------------------------------------------------------------------------------
## Change intercept in the first class (for DM all 5 classes are modeled)
delta <- 2
m1.dm <- m0.dm
m1.dm@coefficients[1, ] <- m0.dm@coefficients[1, ] +
                           c(-delta, rep(delta/4, 4))

alpha0 <- exp(X[phase2,] %*% m0.dm@coefficients)
alpha1 <- exp(X[phase2,] %*% m1.dm@coefficients)

dfun <- function(y, size, mu, log = FALSE) {
  dLog <- ddirmn(t(y), t(mu))
  if (log) dLog else exp(dLog)
}

h <- 2
control <- list(range = seq(nrow(rotaBB))[phase2], h = h,
                pi0 = t(alpha0), pi1 = t(alpha1),
                ret = "value", dfun = dfun)
surv.dm <- categoricalCUSUM(rotaBB, control = control)

## ----echo=FALSE,eval=FALSE---------------------------------------------------------
# matplot(alpha0/rowSums(alpha0),type="l",lwd=3,lty=1,ylim=c(0,1))
# matlines(alpha1/rowSums(alpha1),type="l",lwd=1,lty=2)

## ----ctPlot-simple, echo = FALSE, fig.keep = 'none'--------------------------------
par(mfrow = c(1,2))
surv@multinomialTS <- surv.dm@multinomialTS <- FALSE  # trick plot method ...
plot(surv[,1], col=c(NA,NA,4), ylab = expression(C[t]), ylim = c(0,33),
     xaxis.tickFreq=list("%Y"=atChange, "%m"=atChange),
     xaxis.labelFreq=list("%Y"=atMedian), xaxis.labelFormat="%Y")
abline(h=h, lwd=2, col="darkgrey")
plot(surv.dm[,1], col=c(NA,NA,4), ylab = expression(C[t]), ylim = c(0,33),
     xaxis.tickFreq=list("%Y"=atChange, "%m"=atChange),
     xaxis.labelFreq=list("%Y"=atMedian), xaxis.labelFormat="%Y")
abline(h=h, lwd=2, col="darkgrey")

## ----------------------------------------------------------------------------------
today <- which(epoch(salmNewport) == as.Date("2013-12-23"))
rangeAnalysis <- (today - 4):today
in2013 <- which(isoWeekYear(epoch(salmNewport))$ISOYear == 2013)

algoParameters <- list(range = rangeAnalysis, noPeriods = 10,
                       populationBool = FALSE,
                       b = 4, w = 3, weightsThreshold = 2.58,
                       pastWeeksNotIncluded = 26, pThresholdTrend = 1,
                       thresholdMethod = "nbPlugin", alpha = 0.05,
                       limit54 = c(0, 50))

results <- farringtonFlexible(salmNewport[, c("Baden.Wuerttemberg",
                                             "North.Rhine.Westphalia")],
                              control = algoParameters)

## ----results='asis'----------------------------------------------------------------
start <- isoWeekYear(epoch(salmNewport)[min(rangeAnalysis)])
end <- isoWeekYear(epoch(salmNewport)[max(rangeAnalysis)])
caption <- paste0("Results of the analysis of reported S. Newport ",
                  "counts in two German federal states for the weeks ",
                  start$ISOYear, "-W", start$ISOWeek, " to ",
                  end$ISOYear, "-W", end$ISOWeek,
                  ". Bold red counts indicate weeks with alarms.")
toLatex(results, caption = caption, label = "tableResults",
        ubColumnLabel = "Threshold", include.rownames = FALSE,
        sanitize.text.function = identity)

