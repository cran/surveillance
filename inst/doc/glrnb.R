### R code from vignette source 'glrnb.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: glrnb.Rnw:58-69
###################################################
library(surveillance)
library(xtable)
options(SweaveHooks=list(fig=function() par(mar=c(5,4,4,0),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)))
options(width=70)
set.seed(247)

#####################################################################
# create directory figs if it does not exist
#####################################################################
if(!file.exists("figs/")) dir.create("figs/")



###################################################
### code chunk number 2: glrnb.Rnw:106-108
###################################################
getOption("SweaveHooks")[["fig"]]()
data(shadar)
plot(shadar,main="Number of salmonella hadar cases in Germany 2001-2006")


###################################################
### code chunk number 3: glrnb.Rnw:115-117
###################################################
# Simulate data
simData <- sim.pointSource(length=300,K=0.5,r=0.6,p=0.95)


###################################################
### code chunk number 4: glrnb.Rnw:120-121
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(simData)


###################################################
### code chunk number 5: glrnb.Rnw:154-156
###################################################
getOption("SweaveHooks")[["fig"]]()
survObj <- algo.glrnb(shadar,control=list(range=105:295,alpha=0))
plot(survObj,startyear=2003)


###################################################
### code chunk number 6: glrnb.Rnw:175-178 (eval = FALSE)
###################################################
## control=list(range=range,c.ARL=5, 
##   mu0=NULL, alpha=0, Mtilde=1, M=-1, change="intercept",theta=NULL,
##   dir=c("inc","dec"),ret=c("cases","value"))


###################################################
### code chunk number 7: glrnb.Rnw:187-189 (eval = FALSE)
###################################################
## control=list(range=105:length(shadar$observed))
## algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 8: glrnb.Rnw:195-197 (eval = FALSE)
###################################################
## control=list(range=105:295,alpha=3)
## algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 9: glrnb.Rnw:205-208
###################################################
control=list(range=105:295,alpha=NULL)
surv <- algo.glrnb(shadar,control=control)
surv$control$alpha


###################################################
### code chunk number 10: glrnb.Rnw:219-221 (eval = FALSE)
###################################################
## control=list(range=105:295,mu0=list(S=2,trend=FALSE))
## algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 11: glrnb.Rnw:224-226
###################################################
control=list(range=105:295,mu0=list(S=2,trend=F,refit=T))
surv <- algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 12: glrnb.Rnw:231-233
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(shadar)
with(surv$control,lines(mu0~range,lty=2,lwd=4,col=4))


###################################################
### code chunk number 13: glrnb.Rnw:239-240 (eval = FALSE)
###################################################
## surv$control$mu0Model


###################################################
### code chunk number 14: glrnb.Rnw:245-246
###################################################
estimateGLRNbHook


###################################################
### code chunk number 15: glrnb.Rnw:285-286
###################################################
coef(surv$control$mu0Model$fitted[[1]])


###################################################
### code chunk number 16: glrnb.Rnw:294-297
###################################################
control=list(range=105:295,alpha=0)
surv <- algo.glrnb(disProgObj=shadar,control=control)
table(surv$alarm)


###################################################
### code chunk number 17: glrnb.Rnw:302-306
###################################################
num <- rep(NA)
for (i in 1:6){
num[i] <- table(algo.glrnb(disProgObj=shadar,control=c(control,c.ARL=i))$alarm)[2]
}


###################################################
### code chunk number 18: glrnb.Rnw:331-333 (eval = FALSE)
###################################################
## control=list(range=105:295,theta=0.4)
## algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 19: glrnb.Rnw:338-340 (eval = FALSE)
###################################################
## control=list(range=105:295,theta=NULL)
## algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 20: glrnb.Rnw:348-350
###################################################
control=list(range=105:295,ret="cases",alpha=0)
surv2 <- algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 21: glrnb.Rnw:353-354
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(surv2,startyear=2003)


###################################################
### code chunk number 22: glrnb.Rnw:364-366
###################################################
control=list(range=105:295,ret="cases",dir="dec",alpha=0)
surv3 <- algo.glrnb(disProgObj=shadar,control=control)


###################################################
### code chunk number 23: glrnb.Rnw:369-370
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(surv3,startyear=2003)


