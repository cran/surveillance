### R code from vignette source 'surveillance.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
library("surveillance")
options(SweaveHooks=list(fig=function() par(mar=c(4,4,2,0)+.5)))
options(width=70)

## create directory for plots
dir.create("plots", showWarnings=FALSE)

######################################################################
#Do we need to compute or can we just fetch results
######################################################################
CACHEFILE <- "surveillance-cache.RData"
compute <- !file.exists(CACHEFILE)
message("Doing computations: ", compute)
if(!compute) load(CACHEFILE)


###################################################
### code chunk number 2: surveillance.Rnw:155-157
###################################################
getOption("SweaveHooks")[["fig"]]()
data(k1)
plot(k1,main="Kryptosporidiosis in BW 2001-2005")


###################################################
### code chunk number 3: surveillance.Rnw:217-221
###################################################
set.seed(1234)
sts <- sim.pointSource(p = 0.99, r = 0.5, length = 400,
                       A = 1, alpha = 1, beta = 0, phi = 0,
                       frequency = 1, state = NULL, K = 1.7)


###################################################
### code chunk number 4: surveillance.Rnw:223-224
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(sts)


###################################################
### code chunk number 5: surveillance.Rnw:317-320
###################################################
getOption("SweaveHooks")[["fig"]]()
k1.b660 <- algo.bayes(k1,
  control = list(range = 27:192, b = 0, w = 6, alpha = 0.01))
plot(k1.b660, disease = "k1", firstweek = 1, startyear = 2001)


###################################################
### code chunk number 6: CDC (eval = FALSE)
###################################################
## cntrl <- list(range=300:400,m=1,w=3,b=5,alpha=0.01)
## sts.cdc  <- algo.cdc(sts, control = cntrl)
## sts.farrington <- algo.farrington(sts, control = cntrl)


###################################################
### code chunk number 7: surveillance.Rnw:348-351
###################################################
if (compute) {
cntrl <- list(range=300:400,m=1,w=3,b=5,alpha=0.01)
sts.cdc  <- algo.cdc(sts, control = cntrl)
sts.farrington <- algo.farrington(sts, control = cntrl)
}


###################################################
### code chunk number 8: surveillance.Rnw:354-357
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mfcol=c(1,2))
plot(sts.cdc, legend.opts=NULL)
plot(sts.farrington, legend.opts=NULL)


###################################################
### code chunk number 9: surveillance.Rnw:375-376
###################################################
print(algo.quality(k1.b660))


###################################################
### code chunk number 10: CONTROL
###################################################
control <- list(
  list(funcName = "rki1"), list(funcName = "rki2"),
  list(funcName = "rki3"), list(funcName = "bayes1"),
  list(funcName = "bayes2"), list(funcName = "bayes3"),
  list(funcName = "cdc", alpha=0.05),
  list(funcName = "farrington", alpha=0.05)
)
control <- lapply(control, function(ctrl) {
  ctrl$range <- 300:400; return(ctrl)
})


###################################################
### code chunk number 11: surveillance.Rnw:416-417 (eval = FALSE)
###################################################
## algo.compare(algo.call(sts, control = control))


###################################################
### code chunk number 12: surveillance.Rnw:419-423
###################################################
if (compute) {
  acall <- algo.call(sts, control = control)
}
print(algo.compare(acall), digits = 3)


###################################################
### code chunk number 13: surveillance.Rnw:432-437
###################################################
#Create 10 series
ten <- lapply(1:10,function(x) {
  sim.pointSource(p = 0.975, r = 0.5, length = 400,
                  A = 1, alpha = 1, beta = 0, phi = 0,
                  frequency = 1, state = NULL, K = 1.7)})


###################################################
### code chunk number 14: TENSURV (eval = FALSE)
###################################################
## #Do surveillance on all 10, get results as list
## ten.surv <- lapply(ten,function(ts) {
##   algo.compare(algo.call(ts,control=control))
## })


###################################################
### code chunk number 15: surveillance.Rnw:445-448
###################################################
if (compute) {
#Do surveillance on all 10, get results as list
ten.surv <- lapply(ten,function(ts) {
  algo.compare(algo.call(ts,control=control))
})
}


###################################################
### code chunk number 16: surveillance.Rnw:450-452 (eval = FALSE)
###################################################
## #Average results
## algo.summary(ten.surv)


###################################################
### code chunk number 17: surveillance.Rnw:454-455
###################################################
print(algo.summary(ten.surv), digits = 3)


###################################################
### code chunk number 18: surveillance.Rnw:467-495
###################################################
#Update range in each - cyclic continuation
range = (2*4*52) +  1:length(k1$observed)
control <- lapply(control,function(cntrl) {
  cntrl$range=range;return(cntrl)})

#Auxiliary function to enlarge data
enlargeData <- function(disProgObj, range = 1:156, times = 1){
  disProgObj$observed <- c(rep(disProgObj$observed[range], times),
                           disProgObj$observed)
  disProgObj$state <- c(rep(disProgObj$state[range], times),
                        disProgObj$state)
  return(disProgObj)
}

#Outbreaks
outbrks <- c("m1", "m2", "m3", "m4", "m5", "q1_nrwh", "q2",
             "s1", "s2", "s3", "k1", "n1", "n2", "h1_nrwrp")

#Load and enlarge data.
outbrks <- lapply(outbrks,function(name) {
  data(list=name)
  enlargeData(get(name),range=1:(4*52),times=2)
})

#Apply function to one
one.survstat.surv <- function(outbrk) {
  algo.compare(algo.call(outbrk,control=control))
}


###################################################
### code chunk number 19: surveillance.Rnw:497-498 (eval = FALSE)
###################################################
## algo.summary(lapply(outbrks,one.survstat.surv))


###################################################
### code chunk number 20: surveillance.Rnw:500-504
###################################################
if (compute) {
  res.survstat <- algo.summary(lapply(outbrks,one.survstat.surv))
}
print(res.survstat, digits=3)


###################################################
### code chunk number 21: surveillance.Rnw:514-520
###################################################
if (compute) { # save computed results
    save(list=c("sts.cdc","sts.farrington","acall","res.survstat",
                "ten.surv"),
         file=CACHEFILE)
    tools::resaveRdaFiles(CACHEFILE)
}


