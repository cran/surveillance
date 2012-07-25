
### Define some functions to be S3 generic

animate <- function (object, ...) UseMethod("animate")
R0 <- function (object, ...) UseMethod("R0")
as.epidata <- function (data, ...) UseMethod("as.epidata")
intensityplot <- function (x, ...) UseMethod("intensityplot")

## internal function with methods for "twinSIR" and "simEpidata"
getModel <- function (object, ...) UseMethod("getModel")

## (rather internal) generic with methods for "matrix" and "Spatial"
multiplicity <- function (x, ...) UseMethod("multiplicity")


### Define some function to be S4 generic

if(!isGeneric("plot")) setGeneric("plot", useAsDefault=plot)
if(!isGeneric("aggregate")) setGeneric("aggregate", useAsDefault=aggregate)

## Register "owin" as class in S4 so we can define methods for it
## Note: package "maptools" also registers "owin" as a virtual S4 class by setClass("owin")
if (!isClass("owin")) {
    setOldClass("owin")
}



######################################################################
#Access and replace functions for the "sts" class
######################################################################
#epoch slot
if(!isGeneric("epoch")) setGeneric("epoch", function(x, as.Date=x@epochAsDate) standardGeneric("epoch"))
setMethod("epoch", "sts", function(x, as.Date=x@epochAsDate) {
  if (!as.Date) {
    return(x@epoch)
  } else {
    return(as.Date(x@epoch, origin="1970-01-01"))
  }
})
setGeneric("epoch<-", function(x, value) standardGeneric("epoch<-"))
setReplaceMethod("epoch", "sts", function(x, value) {
 x@epoch <- value
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
###multinomial Time series slot
##control slot
if(!isGeneric("multinomialTS")) setGeneric("multinomialTS", function(x) standardGeneric("multinomialTS"))
setMethod("multinomialTS", "sts", function(x) {
  return(x@multinomialTS)
})
setGeneric("multinomialTS<-", function(x, value) standardGeneric("multinomialTS<-"))
setReplaceMethod("multinomialTS", "sts", function(x, value) {
 x@multinomialTS <- value
 x
})

### neighbourhood matrix slot 
if(!isGeneric("neighbourhood")) setGeneric("neighbourhood", function(x) standardGeneric("neighbourhood"))
setMethod("neighbourhood", "sts", function(x) {
  return(x@neighbourhood)
})
setGeneric("neighbourhood<-", function(x, value) standardGeneric("neighbourhood<-"))
setReplaceMethod("neighbourhood", "sts", function(x, value) {
 x@neighbourhood <- value
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


