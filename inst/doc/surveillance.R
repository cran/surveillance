### R code from vignette source 'surveillance.Rnw'

###################################################
### code chunk number 1: surveillance.Rnw:54-74
###################################################
library(surveillance)
library(xtable)
options(SweaveHooks=list(fig=function() par(mar=c(5,4,4,0),cex.axis=1.5,cex.lab=1.5,cex.main=1.5)))
options(width=70)
set.seed(1234)

#####################################################################
# create directory figs if it does not exist
#####################################################################
if(!file.exists("figs/")) dir.create("figs/")

######################################################################
#Do we need to compute or can we just fetch results
######################################################################
CACHEFILE <- "surveillance-cache.RData"
compute <- !file.exists(CACHEFILE)
#load computed results
if(!compute) load(CACHEFILE)
print(paste("Doing computations:", compute))



###################################################
### code chunk number 2: surveillance.Rnw:169-171
###################################################
getOption("SweaveHooks")[["fig"]]()
data(k1)
plot(k1,main="Kryptosporidiosis in BW 2001-2005")


###################################################
### code chunk number 3: surveillance.Rnw:264-268
###################################################
getOption("SweaveHooks")[["fig"]]()
sts <- sim.pointSource(p = 0.99, r = 0.5, length = 400,
                       A = 1, alpha = 1, beta = 0, phi = 0,
                       frequency = 1, state = NULL, K = 1.7)
plot(sts)


###################################################
### code chunk number 4: surveillance.Rnw:361-364
###################################################
getOption("SweaveHooks")[["fig"]]()
k1.b660 <- algo.bayes(k1, control = list(range = 27:192,b=0,w=6,alpha=0.01))
plot(k1.b660, disease="k1", firstweek = 1, startyear = 2001)



###################################################
### code chunk number 5: CDC (eval = FALSE)
###################################################
## cntrl <- list(range=300:400,m=1,w=3,b=5,alpha=0.01)
## sts.cdc  <- algo.cdc(sts, control = cntrl)
## sts.farrington <- algo.farrington(sts, control = cntrl)


###################################################
### code chunk number 6: surveillance.Rnw:392-395
###################################################
if (compute) {
cntrl <- list(range=300:400,m=1,w=3,b=5,alpha=0.01)
sts.cdc  <- algo.cdc(sts, control = cntrl)
sts.farrington <- algo.farrington(sts, control = cntrl)
}


###################################################
### code chunk number 7: surveillance.Rnw:398-401
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfcol=c(1,2))
plot(sts.cdc, legend.opts=NULL)
plot(sts.farrington, legend.opts=NULL)


###################################################
### code chunk number 8: surveillance.Rnw:419-421
###################################################
print(algo.quality(k1.b660))



###################################################
### code chunk number 9: CONTROL
###################################################
control = list(
  list(funcName = "rki1"),
  list(funcName = "rki2"),
  list(funcName = "rki3"),
  list(funcName = "bayes1"),
  list(funcName = "bayes2"),
  list(funcName = "bayes3"),
  list(funcName = "cdc",alpha=0.05),
  list(funcName = "farrington",alpha=0.05))
control <- lapply(control,function(ctrl) {
  ctrl$range <- 300:400;return(ctrl)})


###################################################
### code chunk number 10: surveillance.Rnw:462-463 (eval = FALSE)
###################################################
## algo.compare(algo.call(sts,control=control))


###################################################
### code chunk number 11: surveillance.Rnw:465-470
###################################################
if (compute) {
  acall <- algo.call(sts,control=control)
}
algo.compare(acall)



###################################################
### code chunk number 12: surveillance.Rnw:482-487
###################################################
#Create 10 series
ten <- lapply(1:10,function(x) {
  sim.pointSource(p = 0.975, r = 0.5, length = 400,
                  A = 1, alpha = 1, beta = 0, phi = 0,
                  frequency = 1, state = NULL, K = 1.7)})


###################################################
### code chunk number 13: TENSURV (eval = FALSE)
###################################################
## #Do surveillance on all 10, get results as list
## ten.surv <- lapply(ten,function(ts) { 
##   algo.compare(algo.call(ts,control=control))
## })


###################################################
### code chunk number 14: surveillance.Rnw:495-498
###################################################
if (compute) {
#Do surveillance on all 10, get results as list
ten.surv <- lapply(ten,function(ts) { 
  algo.compare(algo.call(ts,control=control))
})
}


###################################################
### code chunk number 15: surveillance.Rnw:500-502 (eval = FALSE)
###################################################
## #Average results
## algo.summary(ten.surv)


###################################################
### code chunk number 16: surveillance.Rnw:504-507
###################################################
res <- algo.summary(ten.surv)
res[,5:8] <- round(res[,5:8]*100)/100
res


###################################################
### code chunk number 17: surveillance.Rnw:520-540
###################################################
#Update range in each - cyclic continuation
range = (2*4*52) +  1:length(k1$observed)
control <- lapply(control,function(cntrl) { 
  cntrl$range=range;return(cntrl)})

#Outbreaks
outbrks <- c("m1", "m2", "m3", "m4", "m5", "q1_nrwh", "q2", 
              "s1", "s2", "s3", "k1", "n1", "n2", "h1_nrwrp")

#Load and enlarge data.
outbrks <- lapply(outbrks,function(name) {
  #Load with data
  eval(substitute(data(name),list(name=name)))
  enlargeData(get(name),range=1:(4*52),times=2)
})

#Apply function to one
one.survstat.surv <- function(outbrk) {
  algo.compare(algo.call(outbrk,control=control))
}


###################################################
### code chunk number 18: surveillance.Rnw:542-543 (eval = FALSE)
###################################################
## algo.summary(lapply(outbrks,one.survstat.surv))


###################################################
### code chunk number 19: surveillance.Rnw:545-550
###################################################
if (compute) {
  res.survstat <- algo.summary(lapply(outbrks,one.survstat.surv))
}
print(res.survstat,digits=3)



###################################################
### code chunk number 20: surveillance.Rnw:589-592
###################################################
getOption("SweaveHooks")[["fig"]]()
data(measles.weser)
plot(measles.weser, title="measles in Weser-Ems 2001-2002", 
     xaxis.years=TRUE, startyear= 2001, firstweek=1)


###################################################
### code chunk number 21: surveillance.Rnw:600-601
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(measles.weser,as.one=FALSE,xaxis.years=FALSE)


###################################################
### code chunk number 22: surveillance.Rnw:619-636
###################################################
#####################################################
# measles
if(compute){
  # algo.hhh
  cntrl <- list(linear=TRUE, nseason=1, neighbours=FALSE, 
                negbin="single", lambda=TRUE)
  measles.hhh.neF <- algo.hhh(measles.weser,control=cntrl)

  cntrl$neighbours <- TRUE
  measles.hhh <-algo.hhh(measles.weser,cntrl)

  # algo.hhh.grid
  grid <- create.grid(measles.weser, cntrl, params=list(endemic=c(lower=-0.5, upper=0.5, length=3), epidemic=c(0.1,0.9,5),  negbin=c(0.3,12,5)))
  cat("running a grid search for up to 900 seconds.\n")
  measles.hhh.grid <- algo.hhh.grid(measles.weser, control=cntrl, 
                                    thetastartMatrix=grid, maxTime=900)
}


###################################################
### code chunk number 23: surveillance.Rnw:651-654 (eval = FALSE)
###################################################
## cntrl <- list(linear=TRUE, nseason=1, neighbours=TRUE, 
##               negbin="single", lambda=TRUE)
## measles.hhh <- algo.hhh(measles.weser,control=cntrl)


###################################################
### code chunk number 24: surveillance.Rnw:673-678 (eval = FALSE)
###################################################
## grid <- create.grid(measles.weser, cntrl, 
##    params = list(endemic = c(lower = -0.5, upper = 0.5, length = 3), 
##      epidemic = c(0.1, 0.9, 5),  negbin = c(0.3,12,5)))
## algo.hhh.grid(measles.weser, control=cntrl, thetastartMatrix=grid,
##               maxTime=900)


###################################################
### code chunk number 25: surveillance.Rnw:680-681
###################################################
print(measles.hhh.grid,3)


###################################################
### code chunk number 26: surveillance.Rnw:684-690
###################################################
if (compute) { # save computed results
    save(list=c("sts.cdc","sts.farrington","acall","res","res.survstat",
                "control","ten.surv",ls(pattern="measles")),
         file=CACHEFILE)
    tools::resaveRdaFiles(CACHEFILE)
}


