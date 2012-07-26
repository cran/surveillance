######################################################################
# Init function for stsBP objects. More or less call-through
# to init of sts objects
######################################################################

init.stsBP <- function(.Object, epoch, start=c(2000,1), freq=52, observed, state=0*observed, map=NULL, neighbourhood=NULL, populationFrac=NULL,alarm=NULL,upperbound=NULL, control=NULL,epochAsDate=FALSE,multinomialTS=FALSE,ci=array(NA,dim=c(nrow(observed),ncol(observed),2)),lambda=NULL) {
  
  .Object <- init.sts(.Object, epoch, start, freq, observed, state, map, neighbourhood, populationFrac,alarm,upperbound, control,epochAsDate,multinomialTS)

  #Check that CI matches
  dim.ci <- c(dim(observed),2)
  if (all(dim(ci) == dim.ci)) {
    .Object@ci <- ci
  } else {
    stop("Dimension of confidence interval (ci) (",paste(dim.ci,collapse=","),") is wrong.\n")
  }
  #Assign lambda slot (bootstrap replicates)
  .Object@lambda <- lambda
  
  return(.Object)
}

######################################################################
# Coerce method
######################################################################

setAs("sts", "stsBP", function(from) {
  stsBP <- new("stsBP",
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
               ci=array(NA,dim=c(nrow(from@observed),ncol(from@observed),2)),
               lambda=array(NA,dim=c(nrow(from@observed),ncol(from@observed),2)))
  return(stsBP)
  })


######################################################################
#Methods
######################################################################
setMethod("initialize", "stsBP", init.stsBP)

