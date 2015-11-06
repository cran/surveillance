######################################################################
# Init function for stsNC objects. More or less call-through
# to init of sts objects
######################################################################

init.stsNC <- function(.Object, epoch, start=c(2000,1), freq=52, observed, state=0*observed, map=NULL, neighbourhood=NULL, populationFrac=NULL,alarm=NULL,upperbound=NULL, control=NULL,epochAsDate=FALSE,multinomialTS=FALSE,reportingTriangle=NULL,predPMF=list(),pi=array(NA,dim=c(nrow(observed),ncol(observed),2)),truth=NULL,delayCDF=list(),SR=NULL) {

  #Make the base sts object
  .Object <- init.sts(.Object, epoch, start, freq, observed, state, map, neighbourhood, populationFrac,alarm,upperbound, control,epochAsDate,multinomialTS)

  #Check that CI matches
  dim.pi <- c(dim(observed),2)
  if (all(dim(pi) == dim.pi)) {
    .Object@pi <- pi
  } else {
    stop("Dimension of confidence interval (pi) (",paste(dim.pi,collapse=","),") is wrong.\n")
  }
  #Assign lambda slot (bootstrap replicates)
  if (is.null(reportingTriangle)) {
    .Object@reportingTriangle <- matrix(NA,nrow(observed),ncol(observed))
  } else {
    .Object@reportingTriangle <- reportingTriangle
  }
  .Object@predPMF <- predPMF
  if (is.null(truth)) {
    .Object@truth <- init.sts(.Object, epoch, start, freq, observed, state, map, neighbourhood, populationFrac,alarm,upperbound, control,epochAsDate,multinomialTS)
  } else {
    .Object@truth <- truth
  }
  .Object@delayCDF <- delayCDF
  if (is.null(SR)) {
    .Object@SR <- array(NA,dim=c(nrow(.Object@observed),1,4))
  } else {
    .Object@SR <- SR
  }
  
  return(.Object)
}

######################################################################
# Coerce method
######################################################################

setAs("sts", "stsNC", function(from) {
  stsNC <- new("stsNC",
               epoch=from@epoch,
               freq=from@freq, 
               start=from@start,
               observed=from@observed,
               state = from@state,
               alarm=from@alarm,
               upperbound=from@upperbound,
               neighbourhood=from@neighbourhood,
               populationFrac=from@populationFrac,
               map=from@map,
               control=from@control,
               epochAsDate=from@epochAsDate,
               multinomialTS=from@multinomialTS,
               reportingTriangle=NULL,
               predPMF=list(),
               pi=array(NA,dim=c(nrow(from@observed),ncol(from@observed),2)),
               truth=from,
               delayCDF=list(),
               SR=array(NA,dim=c(nrow(from@observed),1,4)))
  return(stsNC)
  })


######################################################################
#Methods
######################################################################
setMethod("initialize", "stsNC", init.stsNC)

######################################################################
# Create own plotting method for the stsNC class, which starts by
# using the inherited method, but with some additional plotting
# put into the .hookFunSpecial function.
#
# Parameters:
#  same as the for the plot method of sts objects.
######################################################################

setMethod(f="plot", signature=signature(x="stsNC", y="missing"),
          function (x, type = observed ~ time | unit, ...) {

              ## environment of hook function will be set to evaluation 
              ## environment of stsplot_time1() and only then be called
              legend.opts <- lty <- lwd <-
                  "accommodate tools:::.check_code_usage_in_package()"
              
              #Hook function specifically for nowcasting objects.
              nowcastPlotHook <- function() {
                  #Define some colors for the plotting as well as some plot symbols
                  color <- surveillance.options("colors")
                  pchList   <- c(nowSymbol=10)
                  
                  #Prolong line of last observation (this should go into the plot function
                  idx <- nrow(x) - which.max(!is.na(rev(upperbound(x)))) + 1
                  #Continue line from plot - use same style as stsplot_time1
                  lines( idx+c(-0.5,0.5), rep(upperbound(x)[idx,],2),col=col[3],lwd=lwd[3],lty=lty[3])

                  #Add the prediction intervals as bars (where not NA). Conf level
                  #is found in x@control$alpha
                  idxt <- which(apply(x@pi[1:nrow(x),1,],1,function(x) all(!is.na(x))))
                  for (i in idxt) {
                      lines( i+c(-0.3,0.3), rep(x@pi[i,,1],2),lty=1,col=color["piBars"])
                      lines( i+c(-0.3,0.3), rep(x@pi[i,,2],2),lty=1,col=color["piBars"])
                      lines( rep(i,each=2), x@pi[i,,],lty=2,col=color["piBars"])
                  }

                  #Extract now date and date range of the plotting
                  startDate <- epoch(x)[1]

                  #Add "now" symbol on x-axis. Plotting now takes possible temporal aggregation into account.
                  #points(x@control$now-startDate+1,0,pch=pchList["nowSymbol"],col=color["nowSymbol"],cex=1.5)
                  points(x@control$timeDelay(startDate,x@control$now)+1,0,pch=pchList["nowSymbol"],col=color["nowSymbol"],cex=1.5)
                  #Add this to the legend
                  if (!is.null(legend.opts)) {
                      legend(x="topright",c("Now"),pch=pchList["nowSymbol"],col=color["nowSymbol"],bg="white")
                  }
                  
                  return(invisible())
              }

              callNextMethod(x=x, type=type, ..., .hookFuncInheritance=nowcastPlotHook)
          })



