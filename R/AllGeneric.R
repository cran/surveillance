# generate generic functions needed

#Make summary a generic function
setGeneric("summary")

#Conversion of some other functions
if(!isGeneric("plot")) setGeneric("plot", useAsDefault=plot)
if(!isGeneric("aggregate")) setGeneric("aggregate", useAsDefault=aggregate)

######################################################################
#Access and replace functions
######################################################################
#epoch slot
if(!isGeneric("epoch")) setGeneric("epoch", function(x, as.Date=x@epochAsDate) standardGeneric("epoch"))
setMethod("epoch", "sts", function(x, as.Date=x@epochAsDate) {
  if (!as.Date) {
    return(x@week)
  } else {
    return(as.Date(x@week, origin="1970-01-01"))
  }
})
setGeneric("epoch<-", function(x, value) standardGeneric("epoch<-"))
setReplaceMethod("epoch", "sts", function(x, value) {
 x@week <- value
 x
})
# observed slot
if(!isGeneric("observed")) setGeneric("observed", function(x) standardGeneric("observed"))
setMethod("observed", "sts", function(x) {
  return(x@observed)
})
setGeneric("observed<-", function(x, value) standardGeneric("observed<-"))
setReplaceMethod("observed", "sts", function(x, value) {
 x@observed <- value
 x
})
# alarms slot
if(!isGeneric("alarms")) setGeneric("alarms", function(x) standardGeneric("alarms"))
setMethod("alarms", "sts", function(x) {
  return(x@alarm)
})
setGeneric("alarms<-", function(x, value) standardGeneric("alarms<-"))
setReplaceMethod("alarms", "sts", function(x, value) {
 x@alarm <- value
 x
})
# upperbound slot
if(!isGeneric("upperbound")) setGeneric("upperbound", function(x) standardGeneric("upperbound"))
setMethod("upperbound", "sts", function(x) {
  return(x@upperbound)
})
setGeneric("upperbound<-", function(x, value) standardGeneric("upperbound<-"))
setReplaceMethod("upperbound", "sts", function(x, value) {
 x@upperbound <- value
 x
})
# population slot (actually its populationFrac)
if(!isGeneric("population")) setGeneric("population", function(x) standardGeneric("population"))
setMethod("population", "sts", function(x) {
  return(x@populationFrac)
})
setGeneric("population<-", function(x, value) standardGeneric("population<-"))
setReplaceMethod("population", "sts", function(x, value) {
 x@populationFrac <- value
 x
})
##control slot
if(!isGeneric("control")) setGeneric("control", function(x) standardGeneric("control"))
setMethod("control", "sts", function(x) {
  return(x@control)
})
setGeneric("control<-", function(x, value) standardGeneric("control<-"))
setReplaceMethod("control", "sts", function(x, value) {
 x@control <- value
 x
})




######################################################################
#Some access functions similar to matrix/dataframe (definition in sts.R)
######################################################################
if(!isGeneric("nrow")) setGeneric("nrow", useAsDefault=nrow)
if(!isGeneric("ncol")) setGeneric("ncol", useAsDefault=ncol)
if(!isGeneric("colnames")) setGeneric("colnames", useAsDefault=colnames)

#New methods 
if(!isGeneric("as.data.frame")) setGeneric("as.data.frame", useAsDefault=as.data.frame)


