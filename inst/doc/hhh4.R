### R code from vignette source 'hhh4.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: setup
###################################################
library("surveillance")  
library("maptools")
gpclibPermit()
library("Matrix")

options(width=70)
set.seed(247)

#####################################################################
# create directory figs if it does not exist
#####################################################################
if(!file.exists("figs/")) dir.create("figs/")


######################################################################
#Do we need to compute or can we just fetch results
######################################################################
compute <- !file.exists("hhh4-cache.RData")
#load computed results
if(!compute) load("hhh4-cache.RData")
print(paste("Doing computations: ",compute,sep=""))


###################################################
### code chunk number 2: loadInfluMen
###################################################
# load data
data("influMen")
# convert to sts class and print basic information about the time series
print(fluMen <- disProg2sts(influMen))


###################################################
### code chunk number 3: getMen
###################################################
meningo <- fluMen[, "meningococcus"]
dim(meningo)


###################################################
### code chunk number 4: plotfluMen
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(fluMen, type = observed ~ time | unit, # type of plot
             same.scale = FALSE,            # unit-specific ylim ?
             col = "grey"                   # color of bars
             )


###################################################
### code chunk number 5: readInFlu (eval = FALSE)
###################################################
## # read in observed number of cases
## flu.counts <- as.matrix(read.table(system.file("extdata/counts_flu_BYBW.txt", 
##                                       package = "surveillance")))
## # remove 'X' in column names                                      
## colnames(flu.counts) <- substring(colnames(flu.counts),first = 2, last = 5)                                      
## # read in adjacency matrix with elements 1 if two regions share a common border
## nhood <- as.matrix(read.table(system.file("extdata/neighbourhood_BYBW.txt",
##                                       package = "surveillance")))
## # visualize adjacency matrix
## image(Matrix(nhood))


###################################################
### code chunk number 6: nhoodByBw
###################################################
getOption("SweaveHooks")[["fig"]]()
# read in observed number of cases
flu.counts <- as.matrix(read.table(system.file("extdata/counts_flu_BYBW.txt", package = "surveillance")))
# remove 'X' in columnnames                                      
colnames(flu.counts) <- substring(colnames(flu.counts),first=2,last=5)                                      
# read in a shapefile of the districts in Bavaria and Baden-Wuerttemberg
map <- readShapePoly(system.file("shapes/districts_BYBW.shp", package = "surveillance"), IDvar = "id")
# read in adjacency matrix with elements 1 if two regions share a common border
nhood <- as.matrix(read.table(system.file("extdata/neighbourhood_BYBW.txt", package = "surveillance")))
print(image(Matrix(nhood)))


###################################################
### code chunk number 7: fluAsSTS
###################################################
# read in a shapefile of the districts in Bavaria and Baden-Wuerttemberg
map <- readShapePoly(system.file("shapes/districts_BYBW.shp",
                        package = "surveillance"), IDvar = "id")
# read in population fractions
p <- matrix(read.table(system.file("extdata/population_2001-12-31_BYBW.txt",
                          package = "surveillance"), header = TRUE)$popFrac, 
            nrow = nrow(flu.counts), ncol= ncol(flu.counts), byrow = TRUE)
# create sts object
flu <- new("sts", epoch = 1:nrow(flu.counts),
                  observed = flu.counts,
                  start = c(2001, 1),
                  freq = 52,
                  neighbourhood = nhood,
                  map = map,
                  population = p
                  )


###################################################
### code chunk number 8: plot-flu-ByBw
###################################################
getOption("SweaveHooks")[["fig"]]()
par(mar=c(0,0,0,0))
plot(flu[year(flu) == 2001, ],    # select year 2001
     type = observed ~ 1 | unit,  # map of counts aggregated over times t
     labels = FALSE               # suppress region labels in map
     )


###################################################
### code chunk number 9: plot-measles
###################################################
getOption("SweaveHooks")[["fig"]]()
data("measlesDE")
# aggregate into successive bi-weekly periods
measles2w <- aggregate(measlesDE, nfreq = 26)

plot(measles2w, type = observed ~ time,   # plot aggregated over all units i
                   main = "Bi-weekly number of measles cases in Germany", 
                   legend.opts = NULL        # suppress default legend
                   )


###################################################
### code chunk number 10: hhh4 (eval = FALSE)
###################################################
## hhh4(sts, control)


###################################################
### code chunk number 11: controlObj (eval = FALSE)
###################################################
## control = list(
##     ar = list(f = ~ -1),       # formula: exp(u'alpha) * y_i,t-1 
##     ne = list(f = ~ -1,        # formula: exp(x'beta) * sum_j {w_ji * y_j,t-1} 
##               weights = NULL   # matrix with weights w_ji 
##                                # [w_ji = neighbourhood(stsObj) as default]
##               ),              
##     end = list(f = ~ 1,        # formula:  exp(z'gamma) * e_it 
##               offset = NULL    # optional offset e_it 
##               ),
##     family = "Poisson",                # Poisson or NegBin model
##     subset = 2:nrow(stsObj),           # subset of observations to be used 
##                                        # in the fitting process
##     optimizer = list(tech = "nlminb"), # details for optimizer 
##     verbose = FALSE,                   # no progress information is printed
##     start = list(fixed = NULL,         # list with initial values for fixed,
##                  random = NULL,        # random, and
##                  sd.corr = NULL        # variance parameters
##                  ),
##     data = data.frame(t = epoch(sts))  # data.frame,
##                                        # or named list with covariates 
##     )
##            


###################################################
### code chunk number 12: fitMeningo1
###################################################
# specify formula object for endemic component
( f_S1 <- addSeason2formula(f = ~ 1, S = 1, period = 52) )
# fit Poisson model
summary(hhh4(meningo, control = list(end = list(f = f_S1), family = "Poisson")))


###################################################
### code chunk number 13: fitMeningo2
###################################################
result1 <- hhh4(meningo, control = list(end = list(f = f_S1), 
                                        family = "NegBin1")) 


###################################################
### code chunk number 14: fitMeningo3
###################################################
m2 <- list(ar = list(f = ~ 1),     # log(lambda) = alpha
           end = list(f = f_S1), 
           family = "NegBin1",
           # use estimates from previous model as initial values
           start = list(fixed = c(log(0.1),      # initial values for alpha, 
                                  coef(result1)) # and remaining parameters
                        )
           )
# fit model
result2 <- hhh4(meningo, control = m2)
# extract ML estimates
round(coef(result2, se = TRUE,      # also return standard errors
                    idx2Exp = 1     # exponentiate 1st param [-> exp(alpha)]
                    ),2)
# get AIC
AIC(result2)


###################################################
### code chunk number 15: hhh4.Rnw:529-530
###################################################
neighbourhood(fluMen) <- matrix(c(0,0,1,0),2,2)


###################################################
### code chunk number 16: fitFluMen
###################################################
# create formula for endemic component
f.end <- addSeason2formula(f = ~ -1 + fe(1, which = c(TRUE, TRUE)), 
                                           # disease-specific intercepts
                           S = c(3, 1),    # S = 3 for flu, S = 1 for men
                           period = 52)
# specify model
m <- list(ar = list(f = ~ -1 + fe(1, which=c(TRUE, TRUE))), 
          ne = list(f = ~ -1 + fe(1, which=c(FALSE, TRUE))), 
          end = list(f = f.end),
          family = "NegBinM"
          )


# fit model
summary(result <- hhh4(fluMen, control = m))


###################################################
### code chunk number 17: plot-fit_men
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(result, i = 2, col = c("orange", "blue", "grey85"), legend = TRUE)


###################################################
### code chunk number 18: ri (eval = FALSE)
###################################################
## f.end <- ~ -1 + ri(type = "iid", corr = "all")


###################################################
### code chunk number 19: computeFluBYBW
###################################################
# weight matrix w_ji = 1/(No. neighbors of j) if j ~ i, and 0 otherwise
wji <- neighbourhood(flu)/rowSums(neighbourhood(flu))

f.end <- addSeason2formula(f = ~ -1 + ri(type = "iid", corr="all") +I((t-208)/100) , S = 3, period = 52)

model.B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1+ ri(type = "iid", corr="all"), weights = wji),
                 end = list(f = f.end, offset = population(flu)),
                 family = "NegBin1", 
                 verbose = 1, 
                 start=list(fixed=c(-0.9,-1.53,0.56,2.45,2.05,0.33,-0.49,0.21,-0.36,0.21,-0.09),
                            sd.corr=c(-0.02,-0.34,0.68))
                 )


if(compute){   
  # this is time-consuming...
  result.B2 <- hhh4(flu, model.B2)
  s.B2 <- summary(result.B2)

  system.time(pred.B2 <- oneStepAhead(result.B2,tp=nrow(flu)-2*52))
  meanSc.B2 <- colMeans(scores(pred.B2)) 
    
  save(s.B2, meanSc.B2, file="hhh4-cache.RData")
}


###################################################
### code chunk number 20: fitFluBYBW (eval = FALSE)
###################################################
## # weight matrix w_ji = 1/(No. neighbors of j) if j ~ i, and 0 otherwise
## wji <- neighbourhood(flu)/rowSums(neighbourhood(flu))
## 
## # endemic component: iid random effects, linear trend, and S=3 seasonal terms
## f.end <- addSeason2formula(f = ~ -1 + ri(type = "iid", corr="all") + 
##                                 I((t-208)/100), 
##                                 S = 3, 
##                                 period = 52)
## 
## model.B2 <- list(ar = list(f = ~ 1),
##                  ne = list(f = ~ -1+ ri(type = "iid", corr="all"), 
##                            weights = wji),
##                  end = list(f = f.end, offset = population(flu)),
##                  family = "NegBin1"
##                  )
## 
## # fit model
## summary(result.B2 <- hhh4(flu, model.B2))


###################################################
### code chunk number 21: hhh4.Rnw:632-633
###################################################
s.B2


###################################################
### code chunk number 22: oneStepAhead (eval = FALSE)
###################################################
## pred.B2 <- oneStepAhead(result.B2, tp = nrow(flu) - 2 * 52)


###################################################
### code chunk number 23: scores (eval = FALSE)
###################################################
## colMeans(scores(pred.B2)[, c("logs", "rps")])


###################################################
### code chunk number 24: hhh4.Rnw:654-655
###################################################
meanSc.B2[ c("logs", "rps")]


###################################################
### code chunk number 25: createVacc
###################################################
data(MMRcoverageDE)
cardVac1 <- MMRcoverageDE[1:16,3:4]

adjustVac <- function(cardVac, p=0.5,nrow=1){
  card <- cardVac[,1]
  vac <- cardVac[,2]
  vacAdj <- vac*card + p*vac*(1-card)
  return(matrix(vacAdj,nrow=nrow, ncol=length(vacAdj), byrow=TRUE))
}
vac0 <- 1-adjustVac(cardVac1,p=0.5,nrow=measles2w@freq*3)
colnames(vac0) <- colnames(measles2w)


###################################################
### code chunk number 26: hhh4.Rnw:687-688
###################################################
vac0[1:2, 1:5]


###################################################
### code chunk number 27: fitMeasles
###################################################
# endemic component: Intercept + S = 1 sine/cosine pair
f.end <- addSeason2formula(f = ~ 1, S = 1, period = 26)
# autoregressive component: Intercept + vaccination coverage information
model.A0 <- list(ar = list(f = ~ 1 + logVac0),
                 end = list(f = f.end, offset = population(measles2w)),
                 data = list(t = epoch(measles2w), logVac0 = log(vac0)))
# fit model
result.A0 <- hhh4(measles2w, model.A0)    
             
# parameter estimates
round(coef(result.A0, 
           se = TRUE,              # also return standard errors
           amplitudeShift = TRUE   # transform sin/cos terms to 
           ), 2)                   # Amplitude/shift formulation


