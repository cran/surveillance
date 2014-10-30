######################################################################
# Everything belonging to the class sts
######################################################################


######################################################################
# Initialize function for the surveillance time series objects
# -- see documentation in RD files
#
# Dimnames are taken from the observed matrix.
######################################################################

#Ensure that all matrix slots have the same dimnames
fix.dimnames <- function(x) {
  #Make sure all arrays have the same dimnames
  dimnames(x@alarm) <- dimnames(x@state) <- dimnames(x@upperbound) <-
    dimnames(x@populationFrac) <- dimnames(x@observed)
  #Special for neighbourhood
  dimnames(x@neighbourhood) <- list(colnames(x@observed),colnames(x@observed))

  return(x)
}

#constructor function
init.sts <- function(.Object, epoch, start=c(2000,1), freq=52, observed, state=0*observed, map=NULL, neighbourhood=NULL, populationFrac=NULL,alarm=NULL,upperbound=NULL, control=NULL,epochAsDate=FALSE,multinomialTS=FALSE) {

  #If used in constructor
  if(nargs() > 1) {
    #Name handling  
    namesObs <-colnames(observed)
    namesState <- colnames(observed)
    #Ensure observed, state are on matrix form
    observed <- as.matrix(observed)
    state <- as.matrix(state)
  
    #check number of columns of observed and state
    nAreas <- ncol(observed)
    nObs <- nrow(observed)
    if(ncol(observed) != ncol(state)){
      #if there is only one state-vector for more than one area, repeat it
      if(ncol(state)==1)
        state <- ts(matrix(rep(state,nAreas),ncol=nAreas,byrow=FALSE),frequency=frequency(observed))
      else{ 
        cat('wrong dimensions of observed and state \n')
      return(NULL)
      }
    }
    
    #check neighbourhood matrix. 
    if(!is.null(neighbourhood) & (any(dim(neighbourhood) != nAreas))) {
      cat('wrong dimensions of neighbourhood matrix \n')
      return(NULL)
    }
    
    #popFrac
    if (nAreas==1 && (!multinomialTS)) {
        populationFrac <- matrix(1,nrow=nObs, ncol=1)
    } else if (is.null(populationFrac)) {
        populationFrac <- matrix(1/nAreas,nrow=nObs,ncol=nAreas)
    } else if (is.vector(populationFrac, mode="numeric") &&
               length(populationFrac) == nAreas) {
        populationFrac <- matrix(populationFrac, nObs, nAreas, byrow=TRUE)
    }
    
    #labels for observed and state
    if(is.null(namesObs)){
      namesObs <- paste("observed", 1:nAreas, sep="")       
      namesState <- paste("state", 1:nAreas, sep="")  
    }
    
    dimnames(observed) <- list(NULL,namesObs)
    dimnames(state) <- list(NULL,namesState)

    #FIXME: If ncol(observed) is huge then the generated matrix
    #might be beyond memory capacities -> use sparse matrices?
    if (is.null(neighbourhood))
      neighbourhood <- matrix(NA,nrow=ncol(observed),ncol=ncol(observed))
      #diag(neighbourhood) <- 0  #FIXME: shouldn't we always define the diag as 0?
    if (is.null(alarm)) 
      alarm      <- matrix(NA,nrow=dim(observed)[1],ncol=dim(observed)[2])
    if (is.null(upperbound))
      upperbound <- matrix(NA,nrow=dim(observed)[1],ncol=dim(observed)[2])

    ##Assign everything else
    .Object@epoch <- epoch
    .Object@epochAsDate <- epochAsDate
    .Object@multinomialTS <- multinomialTS
    
    if (length(start) == 2) {
      .Object@start <- start
    } else {
      stop("start must be a vector of length two denoting (year, epoch/week/month/idx)")
    }
    
    .Object@freq <- freq
    .Object@state <- state
    .Object@observed <- observed
    
    #It is not possible to assign a null argument to the
    #SpatialPolygonsDataFrame slot. 
    if (!is.null(map)) {
      .Object@map <- map
    }
    
    .Object@neighbourhood <- neighbourhood
    .Object@populationFrac <- populationFrac
    .Object@alarm <- alarm
    .Object@upperbound <- upperbound
    
    if (!is.null(control))
      .Object@control <- control
    
    #Make sure all arrays have the same dimnames
    .Object <- fix.dimnames(.Object)
  }
  
  return(.Object)
}


###########################################################################
# Initialization -- two modes possible: full or just disProg, freq and map
###########################################################################

#Full
setMethod("initialize", "sts", init.sts)

#Partial -- use a disProg object as start and convert it.
disProg2sts <- function(disProgObj, map=NULL) {
  #Ensure that epoch slot is not zero
  if (is.null(disProgObj[["epoch",exact=TRUE]])) {
    myweek <- 1:nrow(as.matrix(disProgObj$observed))
  } else {
    myweek <- disProgObj$week
  }
    
  sts <- new("sts", epoch=myweek, start=disProgObj$start, freq=disProgObj$freq, observed=disProgObj$observed, state = disProgObj$state, map=map, neighbourhood=disProgObj$neighbourhood, populationFrac=disProgObj$populationFrac,alarm=disProgObj$alarm,upperbound=disProgObj$upperbound)
  return(sts)
}

#The reverse action
sts2disProg <- function(sts) {
  disProgObj <- create.disProg(week=sts@epoch, start=sts@start, freq=sts@freq,
                               observed=sts@observed, state=sts@state, neighbourhood=sts@neighbourhood,
                               populationFrac=sts@populationFrac, epochAsDate=sts@epochAsDate)
  #For survRes: alarm=sts@alarm, upperbound=sts@upperbound)
  return(disProgObj)
}


###########################################################################
#Method to aggregate over all units, either the time series is aggregated
#so a new sampling frequency of nfreq units per time slot is obtained.
#The other alternative is to aggregate all units.
#
# Note: The function is not 100% consistent with what the generic
#       aggregate does. 
#
# Warning: In case the aggregation is by unit the upperbound slot is set
#          to NA. Furthermore the MAP object is left as.is, but
#          the object cannot be plotted anymore.
#
# Params:
#   by - a string being either "time" or "unit"
#   nfreq - new sampling frequency if by=="time". If "all" then all
#           time instances are summed.
###########################################################################

setMethod("aggregate", signature(x="sts"), function(x,by="time",nfreq="all",...) {
  
 #Action of aggregation for populationFrac depends on the type 
 binaryTS <- sum( x@populationFrac > 1 ) > 1

  #Aggregate time
  if (by == "time") {
    if (nfreq == "all") {
      howmany <- dim(x@observed)[1]
    } else if (nfreq == x@freq) { # nothing to do
      return(x)
    } else { # nfreq != x@freq
      howmany <- x@freq / nfreq
      if (howmany - ceiling(howmany) != 0)
          stop("nfreq has to be a multiple of x@freq.")
    }
    
    n <- dim(x@observed)[1]
    m <- ceiling(n/howmany)
    new <- rep(1:m,each=howmany)[1:n]
    x@freq <- ifelse(nfreq == "all", howmany, nfreq)
    x@epoch <- 1:m
    
    x@observed <- as.matrix(aggregate(x@observed,by=list(new),sum)[,-1])
    x@state <- as.matrix(aggregate(x@state,by=list(new),sum)[,-1])>0
    x@alarm <- as.matrix(aggregate(x@alarm,by=list(new),sum)[,-1])
    x@upperbound <- as.matrix(aggregate(x@upperbound,by=list(new),sum)[,-1])
    x@populationFrac <- as.matrix(aggregate(x@populationFrac,by=list(new),sum)[,-1])

    #the population fractions need to be recomputed if not a binary ts
    if (!binaryTS) {
      sums <- matrix(rep(apply(x@populationFrac,1,sum),times=ncol(x)),ncol=ncol(x))
      x@populationFrac <-x@populationFrac/sums
    }
  }
  if (by == "unit") {
    #Aggregate units
    x@observed <- as.matrix(apply(x@observed, MARGIN=1, sum))
    x@state <- as.matrix(apply(x@state, MARGIN=1, sum))>0
    x@alarm <- as.matrix(apply(x@alarm, MARGIN=1, sum))>0
    #There is no clever way to aggregate the upperbounds
    x@upperbound <- matrix(NA,ncol=ncol(x@alarm),nrow=nrow(x@alarm))
    x@populationFrac <- as.matrix(apply(x@populationFrac, MARGIN=1, sum))#>0
    x@neighbourhood <- matrix(1,nrow=1,ncol=1)
                     # FIXME: |-wouldn't 0 be more appropriate?
  }

  #validObject(x) #just a check

  return(x)
})
  

#####################################################################
# Miscellaneous access methods
####################################################################

setMethod("nrow", "sts", function(x) return(nrow(x@observed)))
setMethod("ncol", "sts", function(x) return(ncol(x@observed)))
setMethod("dim", "sts", function(x) return(dim(x@observed)))
setMethod("colnames", signature=c(x="sts",do.NULL="missing",prefix="missing"), function(x,do.NULL, prefix) return(colnames(x@observed)))
#Extract which observation within year we have
setGeneric("epochInYear", function(x, ...) standardGeneric("epochInYear"));
setMethod("epochInYear", "sts", function(x,...) {
  #Strptime format strings available as:
  #http://www.opengroup.org/onlinepubs/009695399/functions/strptime.html
  if (x@epochAsDate) {
    epochStr <- switch( as.character(x@freq), "12" = "%m","52" =  "%V","365" = "%j")
    return(as.numeric(formatDate(epoch(x),epochStr)))
  } else {
    return( (x@epoch-1 + x@start[2]-1) %% x@freq + 1)
  }
})
#Extract the corresponding year for each observation using
setGeneric("year", function(x, ...) standardGeneric("year"));
setMethod("year", "sts", function(x,...) {
  if (x@epochAsDate) {
    return(as.numeric(formatDate(epoch(x),"%G")))
  } else {
    ((x@epoch-1 + x@start[2]-1) + (x@freq*x@start[1])) %/% x@freq 
  }
})

#####################################################################
#[-method for accessing the observed, alarm, etc. objects
# new param:
#  normPopulationFrac - normalize population frac
#####################################################################

setMethod("[", "sts", function(x, i, j, ..., drop) {
  #default value for i and j
  if(missing(i)) {i <- min(1,nrow(x@observed)):nrow(x@observed)}
  if(missing(j)) {j <- min(1,ncol(x@observed)):ncol(x@observed)}

  x@epoch <- x@epoch[i]
  x@observed <- x@observed[i,j,drop=FALSE]
  x@state <- x@state[i,j,drop=FALSE]
  x@alarm <- x@alarm[i,j,drop=FALSE]

  x@populationFrac <- x@populationFrac[i,j,drop=FALSE]
  #If not binary TS the populationFrac is normed
  binaryTS <- sum( x@populationFrac > 1 ) > 1
  if (!binaryTS) {
    x@populationFrac <- x@populationFrac / apply(x@populationFrac,MARGIN=1,sum)
   }  
  x@upperbound <- x@upperbound[i,j,drop=FALSE]

  #Neighbourhood matrix
  x@neighbourhood <- x@neighbourhood[j,j,drop=FALSE]

  #Fix the corresponding start entry. it can either be a vector of
  #logicals or a specific index. Needs to work in both cases.
  #Note: This code does not work if we have week 53s!
  if (is.logical(i)) {
    i.min <- which.max(i) #first TRUE entry
  } else {
    i.min <- min(i)
  }
  start <- x@start
  new.sampleNo <- start[2] + i.min - 1
  start.year <- start[1] + (new.sampleNo - 1) %/% x@freq 
  start.sampleNo <- (new.sampleNo - 1) %% x@freq + 1
  x@start <- c(start.year,start.sampleNo)
  ## If !epochAsDate and updating "start" (which seems not really to be
  ## necessary), we have to update epoch, too!
  if (!x@epochAsDate) x@epoch <- x@epoch - i.min + 1  # FIXME: correct?
                                                     
  #FIXME: In case there is a map: Also subset map according to region index j.
  # -> This only makes sense if identical(row.names(map), colnames(observed))
  #    is a property of the sts-class already verified in init.sts
  ## if (length(x@map)>0) {
  ##   ## x@map <- x@map[colnames(x@observed),]
  ##   #take them in the order as in the map to ensure x and x[] are identical
  ##   order <- pmatch(row.names(x@map), colnames(x)) #, row.names(x@map))
  ##   order2 <- pmatch(row.names(x@map)[!is.na(order)], row.names(x@map))
  ##   x@map <- x@map[order2,]
  ## }
  
  #Done
  return(x)
})


#########################################################################
## Plot method ... the type argument specifies what type of plot to make
##
## plot as multivariate time series:  type = observed ~ time | unit 
## plot as map object aggregated over time: type = observed ~ 1 | unit
## new map implementation via: type = observed ~ unit
## the specific plot functions are in separate files (stsplot_*.R)
########################################################################

setMethod("plot", signature(x="sts", y="missing"),
          function (x, type = observed ~ time | unit, ...) {
  
  # catch new implementation of time-aggregate map plot
  if (isTRUE(all.equal(observed ~ unit, type)))
      return(stsplot_space(x, ...))

  #Valid formula?
  valid <- lapply(as.list(type[[3]]), function(i)
                  is.na(pmatch(i,c("1","unit","|","time","*","+"))))
  valid <- all(!unlist(valid))
  obsOk <- (type[[2]] == "observed")
  alarmOk <- (type[[2]] == "alarm")
  if (!valid || !(obsOk | alarmOk))
      stop("Not a valid plot type")

  #Parse the formula, i.e. extract components
  map   <- (length(type[[3]])==3) && (type[[3]][[1]] == "|") && (type[[3]][[2]] == "1")
  time  <- pmatch("time",type[[3]]) > 0
  #All-in-one if type=time+unit -> no, use argument "as.one" for stsplot_time
  #as.one <- all(!is.na(pmatch(c("time","unit"),type[[3]] ))) && is.na(pmatch("|",type[[3]]))

  #No unit dimension?
  justTime <- type[[3]] == "time"

  #space-time plots
  if (map) {
    stsplot_spacetime(x, type, ...)
    return(invisible())
  }
  #time plots
  if (time) {
    if (obsOk) {
      #In case observed ~ time, the units are aggregated
      stsplot_time(if(justTime) aggregate(x,by="unit") else x, ...)
      return(invisible())
    }
    if (alarmOk) {
      stsplot_alarm(x, ...)
      return(invisible())
    }
  }
})


###Validity checking
setValidity("sts", function ( object ) {
  retval <- NULL

  #Check matrix dimensions
  if (!all( dim(object@observed) == dim(object@state)))
    retval <- c( retval , " dimension of observed and state has to match")
  if (!all( dim(object@observed) == dim(object@alarm)))
    retval <- c( retval , " dimension of observed and alarm has to match")
  if (!all( dim(object@observed) == dim(object@upperbound)))
    retval <- c( retval , " dimension of observed and upperbound has to match")
  if (!all( dim(object@observed) == dim(object@populationFrac)))
    retval <- c( retval , " dimension of observed and populationFrac has to match")
  if (!all( rep(dim(object@observed)[2] == dim(object@neighbourhood)))) {
    retval <- c( retval , " dimension of neighbourhood has to be ncols(observed) x ncols(observed)")
  }

  #Check colnames
  if (!all( colnames(object@observed) == colnames(object@state)))
    retval <- c( retval , " colnames of observed and state have to match")
  if (!all( colnames(object@observed) == colnames(object@alarm)))
    retval <- c( retval , " colnames of observed and alarm have to match")
  if (!all( colnames(object@observed) == colnames(object@upperbound)))
    retval <- c( retval , " colnames of observed and upperbound have to match")
  if (!all( colnames(object@observed) == colnames(object@populationFrac)))
    retval <- c( retval , " colnames of observed and populationFrac have to match")
  if (!all( colnames(object@observed) == colnames(object@neighbourhood)))
    retval <- c( retval , " colnames of observed and neighbourhood have to match")

  ## if map is not empty, check that all colnames(object@observed)
  ## are in row.names(object@map)
  if (length(object@map)>0) {
    if (!all(colnames(object@observed) %in% row.names(object@map))) {
      retval <- c( retval , " colnames of observed have to be contained in the row.names of the map.")
    }
  }
  
  if(is.null( retval)) return ( TRUE )
  else return ( retval )
})


##Show/Print method
# show
setMethod( "show", "sts", function( object ){
  cat( "-- An object of class sts -- \n" )
  #cat( "length(week):\t", length(object@epoch),"\n" )
  if (!object@epochAsDate) {
    cat( "freq:\t\t", object@freq,"\n" )
  } else {
    epochStr <- switch( as.character(object@freq), "12" = "%m","52" =  "%V","365" = "%j")
    cat( "freq:\t\t", paste(object@freq," with strptime format string ",epochStr,"\n",sep=""))
  }
  if (!object@epochAsDate) {
    cat( "start:\t\t",object@start,"\n" )
  } else {
    cat( "start:\t\t",paste(epoch(object)[1]),"\n" )
  }
  cat( "dim(observed):\t", dim(object@observed), "\n\n")

  n <- 1
  cat("Head of observed:\n")
  print(head(object@observed,n))
  #cat("Head of state:\n")
  #print(head(object@state,n))
  #cat("Head of alarm:\n")
  #print(head(object@alarm,n))

  if (npoly <- length(object@map)) {
      cat("\nmap:\n")
      print(modifyList(summary(object@map), list(data=NULL))) # no data summary
      cat("Features    :", npoly, "\n")
      if (inherits(object@map, "SpatialPolygonsDataFrame"))
          cat("Data slot   :", ncol(object@map), "variables\n")
  }

  cat("\nhead of neighbourhood:\n")
  print( head(object@neighbourhood,n))
} )


  

 
####################
# toLatex method
####################  

toLatex.sts <- function(object, caption = "",label=" ", columnLabels = NULL,
                        subset = NULL, 
                        alarmPrefix = "\\textbf{\\textcolor{red}{",
                        alarmSuffix = "}}", ubColumnLabel = "UB", ...) {
  # Takes a list of sts objects and outputs a LaTeX Table. 
  # Args:
  #   object: A single sts object; must not be NULL or empty.
  #   caption: A caption for the table. Default is the empty string.
  #   label: A label for the table. Default is the empty string.
  #   columnLabels: A list of labels for each column of the resulting table.
  #   subset: A range of values which should be displayed. If Null, then all 
  #           data in the sts objects will be displayed. Else only a subset of 
  #           data. Therefore range needs to be a numerical vector of indexes
  #           from 1 to length(@observed).
  #   alarmPrefix: A latex compatible prefix string wrapped around a table cell 
  #                iff there is an alarm;i.e. alarm = TRUE
  #   alarmSuffix: A latex compatible suffix string wrapped around a table cell 
  #                iff there is an alarm;i.e. alarm[i,j] = TRUE
  #   ubColumnLabel: The label of the upper bound column; default is "UB".
  #   ...: Variable arguments passed to toLatex.xtable
  # Returns:
  #  An object of class Latex
  
  # Error Handling
  isEmpty <- function(o) is.null(o)
  if (isEmpty(object))
    stop("object must not be null or NA.")
  
  if (is.list(object))
    stop("supplying a list of sts has been removed from the api. Sorry.")
  
  if (!isS4(object) || !is(object, "sts"))
    stop("object must be of type sts from the surveillance package.")
  
  if (!is.character(caption))
    stop("caption must be a character.")
  
  if (!isEmpty(labels) && length(labels) != length(object))
    stop("number of labels differ from the number of sts objects.")
    
  # derive default values
  
  tableLabels <- colnames(object@observed)
  if (!is.null(columnLabels) && 
        length(columnLabels) != ncol(object@observed) * 2 + 2) {
    stop("the number of labels must match the number of columns in the 
         resulting table; i.e. 2 * columns of sts + 2.")
  }
  
  tableCaption <- caption
  tableLabel <- label
  
  epochAsDate <- object@epochAsDate
  epochStr <- switch( as.character(object@freq), 
                      "12" = "month",
                      "52" =  "week",
                      "365" = "day")
  noDataPoints <- nrow(object@observed)
      
  if (epochAsDate) {
    vectorOfDates <- as.Date(object@epoch, origin="1970-01-01")
  } else { 
    firstMonday <- as.POSIXlt(
                              paste(object@start[1],
                                    object@start[2],"1",sep=" ")
                             ,format = "%Y %W %u") 
    dateInc <- switch(as.character(object@freq), 
                          "12" = 30, 
                          "365" = 1, 
                          "52" = 7)
    vectorOfDates <- as.Date(firstMonday) + dateInc * (0:(noDataPoints - 1))
  }
  
  yearColumn <- Map(function(d)isoWeekYear(d)$ISOYear, vectorOfDates)
  
  if (object@freq == 12 )
    monthColumn <- Map(function(d) as.POSIXlt(d)$mon, vectorOfDates)
  
  if (object@freq == 52 )
    weekColumn <- Map(function(d)isoWeekYear(d)$ISOWeek, vectorOfDates)
  
  dataTable <- data.frame(unlist(yearColumn))
  colnames(dataTable) <- "year"
  
  if (object@freq == 12 ) {
    dataTable$month <- unlist(monthColumn)  
  }
  
  if (object@freq == 52 ) {
    dataTable$week <- unlist(weekColumn)
  }
  
  if (object@freq == 365 ) {
    dataTable$day <- unlist(vectorOfDates)
    dataTable <- dataTable[c(2)]
  }
  
  noCols <- ncol(dataTable)
  j <- 1 + noCols
  tableLabelsWithUB <- c()
  
  # I know it is imperative - shame on me
  for (k in 1:(ncol(object@observed))) {
    upperbounds <- round(object@upperbound[,k], 2)
    observedValues <- object@observed[,k]
    alarms <- object@alarm[,k]
    ubCol <- c()
    for (l in 1:length(upperbounds)) {
        if (is.na(upperbounds[l])) {
            ubCol <- c(ubCol, NA)
        } else {
            ubCol <- c(ubCol, upperbounds[l])
            if (!is.na(alarms[l]) && alarms[l]) {
                observedValues[l] <- paste0(alarmPrefix, observedValues[l], alarmSuffix)
            }
        }
    }
    dataTable[,(j)] <- observedValues
    dataTable[,(j + 1)] <- ubCol
    tableLabelsWithUB <- c(tableLabelsWithUB, tableLabels[k])
    tableLabelsWithUB <- c(tableLabelsWithUB, ubColumnLabel)
    j <- j + 2
  }
  
  # remove rows which should not be displayed
  if (is.null(subset))
    subset <- 1:nrow(dataTable)
  else if (min(subset) < 1 || max(subset) > nrow(dataTable))
    stop("'subset' must be a subset of 1:nrow(observed), i.e., 1:",
         nrow(dataTable))
  
  dataTable <- dataTable[subset,]
  
  # prepare everything for xtable
  newColNames <- c(colnames(dataTable)[1:noCols], tableLabelsWithUB)
  if (!is.null(columnLabels)) {
    colnames(dataTable) <- columnLabels
  } else {
    colnames(dataTable) <- newColNames
  }
  xDataTable <- xtable(dataTable, label = tableLabel, caption = tableCaption, digits = c(0))
  toLatex(xDataTable, ...) 
}

setMethod("toLatex", "sts", toLatex.sts)
setMethod("toLatex", "list", toLatex.sts)
##<- FIXME: actually need a formal class for lists of sts objects



