### R code from vignette source 'monitoringCounts.Rnw'

###################################################
### code chunk number 1: SETUP
###################################################
options(width=77)
## create directories for plots and cache
dir.create("plots", showWarnings=FALSE)
dir.create("monitoringCounts-cache", showWarnings=FALSE)
## load packages
library('surveillance')
library('gamlss')


###################################################
### code chunk number 2: monitoringCounts.Rnw:146-147
###################################################
data("salmNewport")


###################################################
### code chunk number 3: test-sts2ts2sts (eval = FALSE)
###################################################
## all.equal(observed(salmNewport),
##           observed(as(as(salmNewport, "ts"), "sts")))


###################################################
### code chunk number 4: monitoringCounts.Rnw:153-156
###################################################
stopifnot(
all.equal(observed(salmNewport),
          observed(as(as(salmNewport, "ts"), "sts")))
)


###################################################
### code chunk number 5: monitoringCounts.Rnw:162-168
###################################################
# This code is the one used for the Salmon et al. (2016) JSS article.
# Using this code all examples from the article can be reproduced.
# computeALL is FALSE to avoid the computationally intensive parts
# of the code (use of simulations to find a threshold value for categoricalCUSUM,
# use of the boda function) but one can set it to TRUE to have it run.
computeALL <- FALSE


###################################################
### code chunk number 6: monitoringCounts.Rnw:170-198
###################################################
# Define plot parameters
#Add lines using grid by a hook function. Use NULL to align with tick marks
hookFunc <- function() { grid(NA,NULL,lwd=1) }
cex.text <- 1.7
cex.axis <- cex.text
cex.main <- cex.text
cex.lab <-  cex.text
cex.leg <- cex.text
line.lwd <- 2#1
stsPlotCol <- c("mediumblue","mediumblue","red2")
alarm.symbol <- list(pch=17, col="red2", cex=2,lwd=3)
#Define list with arguments to use with do.call("legend", legOpts)
legOpts <- list(x="topleft",legend=c(expression(U[t])),bty="n",lty=1,lwd=line.lwd,col=alarm.symbol$col,horiz=TRUE,cex=cex.leg)
#How should the par of each plot look?
par.list <- list(mar=c(6,5,5,5),family="Times")
#Do this once
y.max <- 0
plotOpts <- list(col=stsPlotCol,ylim=c(0,y.max),
                 main='',lwd=c(1,line.lwd,line.lwd),
                 dx.upperbound=0, #otherwise the upperbound line is put 0.5 off
                 cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.main,
                 ylab="No. of reports", xlab="Time (weeks)",lty=c(1,1,1),
                 legend.opts=legOpts,alarm.symbol=alarm.symbol,
                 xaxis.tickFreq=list("%V"=atChange,"%m"=atChange,"%G"=atChange),
                 xaxis.labelFreq=list("%Y"=atMedian),
                 xaxis.labelFormat="%Y",
                 par.list=par.list,hookFunc=hookFunc)



###################################################
### code chunk number 7: NewportPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(aggregate(salmNewport,by="unit")@observed,na.rm=TRUE)
plotOpts2 <- modifyList(plotOpts,list(x=salmNewport,legend.opts=NULL,ylim=c(0,y.max),type = observed ~ time),keep.null=TRUE)
plotOpts2$par.list <- list(mar=c(6,5,0,5),family="Times")
plotOpts2$xaxis.tickFreq <- list("%m"=atChange,"%G"=atChange)
do.call("plot",plotOpts2)


###################################################
### code chunk number 8: monitoringCounts.Rnw:222-225
###################################################
plot(salmNewport, type = observed ~ time,
     xaxis.tickFreq = list("%m" = atChange, "%G" = atChange),
     xaxis.labelFreq = list("%Y" = atMedian), xaxis.labelFormat = "%Y")


###################################################
### code chunk number 9: unitPlot1
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(salmNewport[,2]),observed(salmNewport[,3]),na.rm=TRUE)
plotOpts2 <- modifyList(plotOpts,list(x=salmNewport[,2],legend.opts=NULL,ylim=c(0,y.max)),keep.null=TRUE)
plotOpts2$xaxis.tickFreq <- list("%G"=atChange)
do.call("plot",plotOpts2)


###################################################
### code chunk number 10: unitPlot2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotOpts2 <- modifyList(plotOpts,list(x=salmNewport[,3],legend.opts=NULL,ylim=c(0,y.max)),keep.null=TRUE)
plotOpts2$xaxis.tickFreq <- list("%G"=atChange)
do.call("plot",plotOpts2)


###################################################
### code chunk number 11: EARS
###################################################
in2011 <- which(isoWeekYear(epoch(salmNewport))$ISOYear == 2011)
salmNewportGermany <- aggregate(salmNewport, by = "unit")
control <- list(range = in2011, method = "C1", alpha = 0.05)
surv <- earsC(salmNewportGermany, control = control)
plot(surv)


###################################################
### code chunk number 12: EARSPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(surv),upperbound(surv),na.rm=TRUE)
do.call("plot",modifyList(plotOpts,list(x=surv,ylim=c(0,y.max)),keep.null=TRUE))


###################################################
### code chunk number 13: farHead
###################################################
control1 <- list(range = in2011, noPeriods = 1,
                 b = 4, w = 3, weightsThreshold = 1,
                 pastWeeksNotIncluded = 3, pThresholdTrend = 0.05,
                 thresholdMethod = "delta")
control2 <- list(range = in2011, noPeriods = 10,
                 b = 4, w = 3, weightsThreshold = 2.58,
                 pastWeeksNotIncluded = 26, pThresholdTrend = 1,
                 thresholdMethod = "nbPlugin")


###################################################
### code chunk number 14: monitoringCounts.Rnw:376-377
###################################################
control1$limit54 <- control2$limit54 <- c(0,50)  # for the figure


###################################################
### code chunk number 15: oldVsNewprep
###################################################
salm.farrington <- farringtonFlexible(salmNewportGermany, control1)
salm.noufaily <- farringtonFlexible(salmNewportGermany, control2)


###################################################
### code chunk number 16: farPlot1
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(salm.farrington),upperbound(salm.farrington),observed(salm.noufaily),upperbound(salm.noufaily),na.rm=TRUE)
do.call("plot",modifyList(plotOpts,list(x=salm.farrington,ylim=c(0,y.max))))


###################################################
### code chunk number 17: farPlot2
###################################################
getOption("SweaveHooks")[["fig"]]()
do.call("plot",modifyList(plotOpts,list(x=salm.noufaily,ylim=c(0,y.max))))


###################################################
### code chunk number 18: campyDE (eval = FALSE)
###################################################
## # Load data and create \code{sts}-object
## data("campyDE")
## cam.sts <- sts(epoch=campyDE$date,
##                observed=campyDE$case, state=campyDE$state)
## par(las=1)
## # Plot
## y.max <- max(observed(cam.sts),upperbound(cam.sts),na.rm=TRUE)
## plotOpts3 <- modifyList(plotOpts,list(x=cam.sts,ylab="",legend.opts=NULL,ylim=c(0,y.max),type = observed ~ time),keep.null=TRUE)
## plotOpts3$xaxis.tickFreq <- list("%m"=atChange,"%G"=atChange)
## do.call("plot",plotOpts3)
## par(las=0)
## #mtext(side=2,text="No. of reports",
##      # las=0,line=3, cex=cex.text,family="Times")
## 	 par(family="Times")
## text(-20, 2600, "No. of\n reports", pos = 3, xpd = T,cex=cex.text)
## text(510, 2900, "Absolute humidity", pos = 3, xpd = T,cex=cex.text)
## text(510, 2550, expression(paste("[",g/m^3,"]", sep='')), pos = 3, xpd = T,cex=cex.text)
## lines(campyDE$hum*50, col="white", lwd=2)
## axis(side=4, at=seq(0,2500,by=500),labels=seq(0,50,by=10),las=1,cex.lab=cex.text, cex=cex.text,cex.axis=cex.text,pos=length(epoch(cam.sts))+20)
## #mtext(side=4,text=expression(paste("Absolute humidity [ ",g/m^3,"]", sep='')),
##      # las=0,line=1, cex=cex.text,family="Times")
## 


###################################################
### code chunk number 19: campyDE
###################################################
getOption("SweaveHooks")[["fig"]]()
# Load data and create \code{sts}-object
data("campyDE")
cam.sts <- sts(epoch=campyDE$date,
               observed=campyDE$case, state=campyDE$state)
par(las=1)
# Plot
y.max <- max(observed(cam.sts),upperbound(cam.sts),na.rm=TRUE)
plotOpts3 <- modifyList(plotOpts,list(x=cam.sts,ylab="",legend.opts=NULL,ylim=c(0,y.max),type = observed ~ time),keep.null=TRUE)
plotOpts3$xaxis.tickFreq <- list("%m"=atChange,"%G"=atChange)
do.call("plot",plotOpts3)
par(las=0)
#mtext(side=2,text="No. of reports",
     # las=0,line=3, cex=cex.text,family="Times")
	 par(family="Times")
text(-20, 2600, "No. of\n reports", pos = 3, xpd = T,cex=cex.text)
text(510, 2900, "Absolute humidity", pos = 3, xpd = T,cex=cex.text)
text(510, 2550, expression(paste("[",g/m^3,"]", sep='')), pos = 3, xpd = T,cex=cex.text)
lines(campyDE$hum*50, col="white", lwd=2)
axis(side=4, at=seq(0,2500,by=500),labels=seq(0,50,by=10),las=1,cex.lab=cex.text, cex=cex.text,cex.axis=cex.text,pos=length(epoch(cam.sts))+20)
#mtext(side=4,text=expression(paste("Absolute humidity [ ",g/m^3,"]", sep='')),
     # las=0,line=1, cex=cex.text,family="Times")




###################################################
### code chunk number 20: campyDE-simple
###################################################
data("campyDE")
cam.sts <- sts(epoch = campyDE$date,
               observed = campyDE$case, state = campyDE$state)
plot(cam.sts, col = "mediumblue")
lines(campyDE$hum * 50, col = "white", lwd = 2)
axis(4, at = seq(0, 2500, by = 500), labels = seq(0, 50, by = 10))


###################################################
### code chunk number 21: boda (eval = FALSE)
###################################################
## library("INLA")
## rangeBoda <- which(epoch(cam.sts) >= as.Date("2007-01-01"))
## control.boda <- list(range = rangeBoda, X = NULL, trend = TRUE,
##                      season = TRUE, prior = "iid", alpha = 0.025,
##                      mc.munu = 10000, mc.y = 1000,
##                      samplingMethod = "marginals")
## boda <- boda(cam.sts, control = control.boda)


###################################################
### code chunk number 22: boda-cache
###################################################
if (computeALL) {
##hoehle 2018-07-18: changed code to use NICELOOKINGboda, but that's iid. Reason:
##The option 'rw1' currently crashes INLA.
library("INLA")
rangeBoda <- which(epoch(cam.sts) >= as.Date("2007-01-01"))
control.boda <- list(range = rangeBoda, X = NULL, trend = TRUE,
                     season = TRUE, prior = "iid", alpha = 0.025,
                     mc.munu = 10000, mc.y = 1000,
                     samplingMethod = "marginals")
boda <- boda(cam.sts, control = control.boda)
save(list = c("boda", "control.boda", "rangeBoda"),
     file = "monitoringCounts-cache/boda.RData")
} else {
  load("monitoringCounts-cache/boda.RData")
}


###################################################
### code chunk number 23: boda2 (eval = FALSE)
###################################################
## covarNames <- c("l1.hum", "l2.hum", "l3.hum", "l4.hum",
##                 "newyears", "christmas", "O104period")
## control.boda2 <- modifyList(control.boda,
##                             list(X = campyDE[, covarNames], season = FALSE))
## boda.covars <- boda(cam.sts, control = control.boda2)


###################################################
### code chunk number 24: boda2-cache
###################################################
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


###################################################
### code chunk number 25: bPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(boda.covars),upperbound(boda.covars),na.rm=TRUE)
plotOpts2 <- modifyList(plotOpts,list(x=boda.covars,ylim=c(0,y.max)),keep.null=TRUE)
plotOpts2$xaxis.tickFreq <- list("%m"=atChange,"%G"=atChange)
do.call("plot",plotOpts2)


###################################################
### code chunk number 26: boda3
###################################################
control.far <- list(range=rangeBoda,b=4,w=5,alpha=0.025*2)
far <- farrington(cam.sts,control=control.far)
#Both farringtonFlexible and algo.bayes uses a one-sided interval just as boda.
control.far2 <-modifyList(control.far,list(alpha=0.025))
farflex <- farringtonFlexible(cam.sts,control=control.far2)
bayes <- suppressWarnings(bayes(cam.sts,control=control.far2))


###################################################
### code chunk number 27: boda4
###################################################
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


###################################################
### code chunk number 28: alarmplot (eval = FALSE)
###################################################
## # Make an artifical object containing two columns - one with the boda output
## # and one with the farrington output
## 
## cam.surv <- combineSTS(list(boda.covars=boda.covars,boda=boda,bayes=bayes,
##                             farrington=far,farringtonFlexible=farflex))
## par(mar=c(4,8,2.1,2),family="Times")
## plot(cam.surv,type = alarm ~ time,lvl=rep(1,ncol(cam.surv)),
##      alarm.symbol=list(pch=17, col="red2", cex=1,lwd=3),
##      cex.axis=1,xlab="Time (weeks)",cex.lab=1,xaxis.tickFreq=list("%m"=atChange,"%G"=atChange),xaxis.labelFreq=list("%G"=at2ndChange),
##      xaxis.labelFormat="%G")


###################################################
### code chunk number 29: alarmplot
###################################################
getOption("SweaveHooks")[["fig"]]()
# Make an artifical object containing two columns - one with the boda output
# and one with the farrington output

cam.surv <- combineSTS(list(boda.covars=boda.covars,boda=boda,bayes=bayes,
                            farrington=far,farringtonFlexible=farflex))
par(mar=c(4,8,2.1,2),family="Times")
plot(cam.surv,type = alarm ~ time,lvl=rep(1,ncol(cam.surv)),
     alarm.symbol=list(pch=17, col="red2", cex=1,lwd=3),
     cex.axis=1,xlab="Time (weeks)",cex.lab=1,xaxis.tickFreq=list("%m"=atChange,"%G"=atChange),xaxis.labelFreq=list("%G"=at2ndChange),
     xaxis.labelFormat="%G")



###################################################
### code chunk number 30: glrnb
###################################################
phase1 <- which(isoWeekYear(epoch(salmNewportGermany))$ISOYear < 2011)
phase2 <- in2011
control <- list(range = phase2, c.ARL = 4, theta = log(2), ret = "cases",
                mu0 = list(S = 1, trend = TRUE, refit = FALSE))
salmGlrnb <- glrnb(salmNewportGermany, control = control)


###################################################
### code chunk number 31: glrnbPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(salmGlrnb),upperbound(salmGlrnb),na.rm=TRUE)
do.call("plot",modifyList(plotOpts,list(x=salmGlrnb,ylim=c(0,y.max))))


###################################################
### code chunk number 32: cat
###################################################
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


###################################################
### code chunk number 33: catPlot1 (eval = FALSE)
###################################################
## y.max <- max(observed(salmHospitalized)/population(salmHospitalized),upperbound(salmHospitalized)/population(salmHospitalized),na.rm=TRUE)
## plotOpts2 <- modifyList(plotOpts,list(x=salmHospitalized,legend.opts=NULL,ylab="",ylim=c(0,y.max)),keep.null=TRUE)
## plotOpts2$xaxis.tickFreq <- list("%G"=atChange,"%m"=atChange)
## plotOpts2$par.list <- list(mar=c(6,5,5,5),family="Times",las=1)
## do.call("plot",plotOpts2)
## lines(salmHospitalized@populationFrac/4000,col="grey80",lwd=2)
## lines(campyDE$hum*50, col="white", lwd=2)
## axis(side=4, at=seq(0,2000,by=500)/4000,labels=as.character(seq(0,2000,by=500)),las=1, cex=2,cex.axis=1.5,pos=length(observed(salmHospitalized))+20)
## par(family="Times")
## text(-20, 0.6, "Proportion", pos = 3, xpd = T,cex=cex.text)
## text(520, 0.6, "Total number of \n reported cases", pos = 3, xpd = T,cex=cex.text)


###################################################
### code chunk number 34: catPlot1
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(salmHospitalized)/population(salmHospitalized),upperbound(salmHospitalized)/population(salmHospitalized),na.rm=TRUE)
plotOpts2 <- modifyList(plotOpts,list(x=salmHospitalized,legend.opts=NULL,ylab="",ylim=c(0,y.max)),keep.null=TRUE)
plotOpts2$xaxis.tickFreq <- list("%G"=atChange,"%m"=atChange)
plotOpts2$par.list <- list(mar=c(6,5,5,5),family="Times",las=1)
do.call("plot",plotOpts2)
lines(salmHospitalized@populationFrac/4000,col="grey80",lwd=2)
lines(campyDE$hum*50, col="white", lwd=2)
axis(side=4, at=seq(0,2000,by=500)/4000,labels=as.character(seq(0,2000,by=500)),las=1, cex=2,cex.axis=1.5,pos=length(observed(salmHospitalized))+20)
par(family="Times")
text(-20, 0.6, "Proportion", pos = 3, xpd = T,cex=cex.text)
text(520, 0.6, "Total number of \n reported cases", pos = 3, xpd = T,cex=cex.text)


###################################################
### code chunk number 35: catbis
###################################################
vars <- c( "y", "n", "t", "epochInPeriod", "weekNumber")
m.bbin <- gamlss(cbind(y, n-y) ~ 1 + t
                 + sin(2 * pi * epochInPeriod) + cos(2 * pi * epochInPeriod)
                 + sin(4 * pi * epochInPeriod) + cos(4 * pi * epochInPeriod)
                 + I(weekNumber == 1) + I(weekNumber == 2),
                 sigma.formula =~ 1,
                 family = BB(sigma.link = "log"),
                 data = salmHospitalized.df[phase1, vars])


###################################################
### code chunk number 36: cat2
###################################################
R <- 2
h <- 2
pi0 <- predict(m.bbin, newdata = salmHospitalized.df[phase2, vars],
               type = "response")
pi1 <- plogis(qlogis(pi0) + log(R))
pi0m <- rbind(pi0, 1 - pi0)
pi1m <- rbind(pi1, 1 - pi1)


###################################################
### code chunk number 37: cat2bis
###################################################
populationHosp <- unname(cbind(
  population(salmHospitalized),
  population(salmHospitalized)))
observedHosp <- cbind(
  "Yes" = as.vector(observed(salmHospitalized)),
  "No" = as.vector(population(salmHospitalized) - observed(salmHospitalized)))
salmHospitalized.multi <- sts(
  freq = 52, start = c(2004, 1), epoch = epoch(salmHospitalized),
  observed = observedHosp, population = populationHosp,
  multinomialTS = TRUE)


###################################################
### code chunk number 38: cat2ter
###################################################
dBB.cusum <- function(y, mu, sigma, size, log = FALSE) {
  dBB(if (is.matrix(y)) y[1,] else y,
      if (is.matrix(y)) mu[1,] else mu,
      sigma = sigma, bd = size, log = log)
}


###################################################
### code chunk number 39: cat3
###################################################
controlCat <- list(range = phase2, h = 2, pi0 = pi0m, pi1 = pi1m,
                   ret = "cases", dfun = dBB.cusum)
salmHospitalizedCat <- categoricalCUSUM(salmHospitalized.multi,
                                        control = controlCat,
                                        sigma = exp(m.bbin$sigma.coef))


###################################################
### code chunk number 40: monitoringCounts.Rnw:1022-1023
###################################################
h.grid <- seq(1, 10, by = 0.5)


###################################################
### code chunk number 41: cath (eval = FALSE)
###################################################
## simone <- function(sts, h) {
##   y <- rBB(length(phase2), mu = pi0m[1, , drop = FALSE],
##            bd = population(sts)[phase2, ], sigma = exp(m.bbin$sigma.coef))
##   observed(sts)[phase2, ] <- cbind(y, population(sts)[phase2, 1] - y)
##   one.surv <- categoricalCUSUM(sts,
##                                control = modifyList(controlCat, list(h = h)),
##                                sigma = exp(m.bbin$sigma.coef))
##   return(any(alarms(one.surv)[, 1]))
## }
## set.seed(123)
## nSims <- 1000
## pMC <- sapply(h.grid, function(h) {
##   mean(replicate(nSims, simone(salmHospitalized.multi, h)))
## })
## 
## pMarkovChain <- sapply(h.grid, function(h) {
##   TA <- LRCUSUM.runlength(mu = pi0m[1,,drop = FALSE],
##                           mu0 = pi0m[1,,drop = FALSE],
##                           mu1 = pi1m[1,,drop = FALSE],
##                           n = population(salmHospitalized.multi)[phase2, ],
##                           h = h, dfun = dBB.cusum,
##                           sigma = exp(m.bbin$sigma.coef))
##   return(tail(TA$cdf, n = 1))
## })


###################################################
### code chunk number 42: cath-cache
###################################################
if (computeALL) {
simone <- function(sts, h) {
  y <- rBB(length(phase2), mu = pi0m[1, , drop = FALSE],
           bd = population(sts)[phase2, ], sigma = exp(m.bbin$sigma.coef))
  observed(sts)[phase2, ] <- cbind(y, population(sts)[phase2, 1] - y)
  one.surv <- categoricalCUSUM(sts,
                               control = modifyList(controlCat, list(h = h)),
                               sigma = exp(m.bbin$sigma.coef))
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


###################################################
### code chunk number 43: catF
###################################################
getOption("SweaveHooks")[["fig"]]()
y.max <- max(observed(salmHospitalizedCat[,1])/population(salmHospitalizedCat[,1]),upperbound(salmHospitalizedCat[,1])/population(salmHospitalizedCat[,1]),na.rm=TRUE)
plotOpts3 <- modifyList(plotOpts,list(x=salmHospitalizedCat[,1],ylab="Proportion",ylim=c(0,y.max)))
plotOpts3$legend.opts <- list(x="top",bty="n",legend=c(expression(U[t])),lty=1,lwd=line.lwd,col=alarm.symbol$col,horiz=TRUE,cex=cex.leg)
do.call("plot",plotOpts3)


###################################################
### code chunk number 44: catARL
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(6,5,5,5),family="Times")
matplot(h.grid, cbind(pMC,pMarkovChain),type="l",ylab=expression(P(T[A] <= 56 * "|" * tau * "=" * infinity)),xlab="Threshold h",col=1,cex=cex.text,
cex.axis =cex.text,cex.lab=cex.text)
prob <- 0.1
lines(range(h.grid),rep(prob,2),lty=5,lwd=2)
axis(2,at=prob,las=1,cex.axis=0.7,labels=FALSE)
par(family="Times")
legend(4,0.08,c("Monte Carlo","Markov chain"), lty=1:2,col=1,cex=cex.text,bty="n")


###################################################
### code chunk number 45: ROTAPLOT
###################################################
data("rotaBB")
plot(rotaBB)


###################################################
### code chunk number 46: monitoringCounts.Rnw:1116-1124
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(5.1,20.1,4.1,0),family="Times")
plot(rotaBB,xlab="Time (months)",ylab="",
     col="mediumblue",cex=cex.text,cex.lab=cex.text,cex.axis=cex.text,cex.main=cex.text,
     xaxis.tickFreq=list("%G"=atChange),
     xaxis.labelFreq=list("%G"=at2ndChange),
     xaxis.labelFormat="%G")
par(las=0,family="Times")
mtext("Proportion of reported cases", side=2, line=19, cex=1)


###################################################
### code chunk number 47: monitoringCounts.Rnw:1132-1159
###################################################
# Select a palette for drawing
pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")
#= RColorBrewer::brewer.pal("Set1",n=ncol(rotaBB))

# Show time series of monthly proportions (matplot does not work with dates)
plotTS <- function(prop=TRUE) {
  for (i in 1:ncol(rotaBB)) {
    fun <- if (i==1) plot else lines
    if (!prop) {
      fun(epoch(rotaBB),observed(rotaBB)[,i],type="l",xlab="Time (months)",ylab="Reported cases",ylim=c(0,max(observed(rotaBB))),col=pal[i],lwd=2)
    } else {
      fun(epoch(rotaBB),observed(rotaBB)[,i,drop=FALSE]/rowSums(observed(rotaBB)),type="l",xlab="Time (months)",ylab="Proportion of reported cases",ylim=c(0,max(observed(rotaBB)/rowSums(observed(rotaBB)))),col=pal[i],lwd=2)
    }
  }
  # Add legend
  axis(1,at=as.numeric(epoch(rotaBB)),labels=FALSE,tck=-0.01)
  legend(x="left",colnames(rotaBB),col=pal,lty=1,lwd=2,bg="white")
}

# plotTS(prop=TRUE)
# Show absolute cases
plotTS(prop=FALSE)

# Even easier
rotaBB.copy <- rotaBB ; rotaBB.copy@multinomialTS <- FALSE
plot(rotaBB.copy)



###################################################
### code chunk number 48: monitoringCounts.Rnw:1165-1179
###################################################
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


###################################################
### code chunk number 49: monitoringCounts.Rnw:1185-1191
###################################################
m1 <- m0

m1@coefficients[1, ] <- m0@coefficients[1, ] + log(2)

pi0 <- t(predict(m0, newdata = X[phase2, ])[, reorder])
pi1 <- t(predict(m1, newdata = X[phase2, ])[, reorder])


###################################################
### code chunk number 50: CATCUSUM
###################################################
dfun <- function(y, size, mu, log = FALSE) {
  dmultinom(x = y, size = size, prob = mu, log = log)
}

h <- 2  # threshold for the CUSUM statistic
control <- list(range = seq(nrow(rotaBB))[phase2], h = h, pi0 = pi0,
                pi1 = pi1, ret = "value", dfun = dfun)
surv <- categoricalCUSUM(rotaBB,control=control)


###################################################
### code chunk number 51: monitoringCounts.Rnw:1208-1210 (eval = FALSE)
###################################################
## alarmDates <- epoch(surv)[which(alarms(surv)[,1]==1)]
## format(alarmDates,"%b %Y")


###################################################
### code chunk number 52: CATCUSUMMC (eval = FALSE)
###################################################
## #Number of MC samples
## nSamples <- 1e4
## 
## #Do MC
## simone.stop <- function(sts, control) {
##   phase2Times <- seq(nrow(sts))[phase2]
##   #Generate new phase2 data from the fitted in control model
##   y <- sapply(1:length(phase2Times), function(i) {
##     rmultinom(n=1, prob=pi0[,i],size=population(sts)[phase2Times[i],1])
##   })
##   observed(sts)[phase2Times,] <- t(y)
##   one.surv <- categoricalCUSUM(sts, control=control)
##   #compute P(S<=length(phase2))
##   return(any(alarms(one.surv)[,1]>0))
## }
## 
## set.seed(1233)
## rlMN <- replicate(nSamples, simone.stop(rotaBB, control=control))
## mean(rlMN)  # 0.5002


###################################################
### code chunk number 53: monitoringCounts.Rnw:1240-1243
###################################################
m0.dm <- MGLMreg(as.matrix(rotaBB.df[phase1, 1:5]) ~ -1 + X[phase1, ],
                dist = "DM")
c(m0@AIC, m0.dm@AIC)


###################################################
### code chunk number 54: monitoringCounts.Rnw:1250-1269
###################################################
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


###################################################
### code chunk number 55: monitoringCounts.Rnw:1272-1274 (eval = FALSE)
###################################################
## matplot(alpha0/rowSums(alpha0),type="l",lwd=3,lty=1,ylim=c(0,1))
## matlines(alpha1/rowSums(alpha1),type="l",lwd=1,lty=2)


###################################################
### code chunk number 56: ctPlot1
###################################################
getOption("SweaveHooks")[["fig"]]()
surv@observed[,1] <- 0
surv@multinomialTS <- FALSE
surv.dm@observed[,1] <- 0
surv.dm@multinomialTS <- FALSE
y.max <- max(observed(surv.dm[,1]),upperbound(surv.dm[,1]),observed(surv[,1]),upperbound(surv[,1]),na.rm=TRUE)
plotOpts3 <- modifyList(plotOpts,list(x=surv[,1],ylim=c(0,y.max),ylab=expression(C[t]),xlab=""))
plotOpts3$legend.opts <- list(x="topleft",bty="n",legend="R",lty=1,lwd=line.lwd,col=alarm.symbol$col,horiz=TRUE,cex=cex.leg)
do.call("plot",plotOpts3)
lines( c(0,1e99), rep(h,2),lwd=2,col="darkgray",lty=1)
par(family="Times")
mtext(side=1,text="Time (weeks)",
      las=0,line=3, cex=cex.text)


###################################################
### code chunk number 57: ctPlot2
###################################################
getOption("SweaveHooks")[["fig"]]()
plotOpts3 <- modifyList(plotOpts,list(x=surv.dm[,1],ylim=c(0,y.max),ylab=expression(C[t]),xlab=""))
plotOpts3$legend.opts <- list(x="topleft",bty="n",legend="R",lty=1,lwd=line.lwd,col=alarm.symbol$col,horiz=TRUE,cex=cex.text)
y.max <- max(observed(surv.dm[,1]),upperbound(surv.dm[,1]),observed(surv[,1]),upperbound(surv[,1]),na.rm=TRUE)
do.call("plot",plotOpts3)
lines( c(0,1e99), rep(h,2),lwd=2,col="darkgray",lty=1)
par(family="Times")
mtext(side=1,text="Time (weeks)",
      las=0,line=3, cex=cex.text)


###################################################
### code chunk number 58: monitoringCounts.Rnw:1400-1414
###################################################
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


###################################################
### code chunk number 59: monitoringCounts.Rnw:1417-1427
###################################################
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


