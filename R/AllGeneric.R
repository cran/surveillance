# generate generic functions needed

#Make summary a generic function
setGeneric("summary")

#Conversion of some other functions
if(!isGeneric("plot")) setGeneric("plot", useAsDefault=plot)
if(!isGeneric("aggregate")) setGeneric("aggregate", useAsDefault=aggregate)

#Some access functions
if(!isGeneric("nrow")) setGeneric("nrow", useAsDefault=nrow)
if(!isGeneric("ncol")) setGeneric("ncol", useAsDefault=ncol)
if(!isGeneric("colnames")) setGeneric("colnames", useAsDefault=colnames)

