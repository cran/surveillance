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
               multinomialTS=from@epochAsDate,
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

