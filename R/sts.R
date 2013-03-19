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
    if (is.null(populationFrac)) {
      populationFrac <- matrix(1/nAreas,nrow=nObs,ncol=nAreas)
    }
    if (nAreas ==1 & (!multinomialTS)){
      populationFrac <- matrix(1,nrow=nObs, ncol=1)
    }
    
    #labels for observed and state
    if(is.null(namesObs)){
      namesObs <- paste("observed", 1:nAreas, sep="")       
      namesState <- paste("state", 1:nAreas, sep="")  
    }
    
    dimnames(observed) <- list(NULL,namesObs)
    dimnames(state) <- list(NULL,namesState)

    #Problem: If ncol(observed) is huge then the generated matrix
    #might be beyond memory capacities -> FIXME: use sparse matrices
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
    } else {
      if (nfreq != x@freq) {
        howmany <- x@freq / nfreq
        if (howmany - ceiling(howmany) != 0) { stop("Error: nfreq has to be a multiple of x@freq.")}
      }
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

  #Fix the corresponding start entry. i can either be a vector of
  #logicals or specific index. Needs to work in both cases.
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

  #Save time by not allocating a new object
  #res <- new("sts",epoch=week, freq=x@freq, start=start,observed=observed,state=state,alarm=alarm,upperbound=upperbound,neighbourhood=neighbourhood,populationFrac=populationFrac,map=x@map,control=x@control)

  ## FIXME: also subset map according to region index j? Here's the code:
  ##x@map <- x@map[colnames(x@observed),]
  
  return(x)
})

#########################################################################
# Plot method ... the type argument specifies what type of plot
# to make.
#
# plot as multivariate time series:  type = observed ~ time | unit 
# plot as map object aggregated over time: type = observed ~ 1 | unit
########################################################################
setMethod("plot", signature(x="sts", y="missing"), function(x, y, type,...) {
  if (missing(type)) type = observed ~ time | unit
  
  #Parse the formula, i.e. extract components
  obsOk <- (type[[2]] == "observed")
  alarmOk <- (type[[2]] == "alarm")
  map   <- (length(type[[3]])==3) && (type[[3]][[1]] == "|") && (type[[3]][[2]] == "1")
  time  <- pmatch("time",type[[3]]) > 0

  #Valid formula?
  valid <- lapply(as.list(type[[3]]),function(i) is.na(pmatch(i,c("1","unit","|","time","*","+"))))
  valid <- all(!unlist(valid))

  #No unit dimenstion?
  justTime <- type[[3]] == "time"
  
  if (!(obsOk | alarmOk) | !valid) {
    stop("Not a valid plot type.")
  }


  #space-time plots
  if (map) {
    plot.sts.spacetime(x,type,...)
    return(invisible())
  }
  #time plots
  if (time) {
    if (obsOk) {
      #In case observed ~ time, the units are aggregated
      plot.sts.time( if(justTime) aggregate(x,by="unit") else x,type,...)
      return(invisible())
    }
    if (alarmOk) {
      plot.sts.alarm(x,...)
      return(invisible())
    }
  }
})

######################################################################
# Helper function to merge two lists taken from the RCurl package
######################################################################

merge.list <- function (x, y, ...) 
{
    if (length(x) == 0) 
        return(y)
    if (length(y) == 0) 
        return(x)
    i = match(names(y), names(x))
    i = is.na(i)
    if (any(i)) 
        x[names(y)[which(i)]] = y[which(i)]
    return(x)
}

##########################################################################
# Plot functions
#
# colors - c( fill color of polygons, line color of polygons, upperbound)
##########################################################################

addFormattedXAxis <- function(x, epochsAsDate, observed, firstweek,xaxis.units,cex) {
  #Declare commonly used variables.
  startyear <-  x@start[1]
  
  if (x@freq ==52) {
    if (!epochsAsDate) {
      # At which indices to put the "at" tick label. This will
      # be exactly those week numbers where the new quarter begins: 1, 14, 27 and 40 + i*52.
      # Note that week number and index is not the same due to the "firstweek" argument
      weeks <- 1:length(observed) + (firstweek-1)
      noYears <- ceiling(max(weeks)/52)
      quarterStarts <- rep( (0:(noYears))*52, each=4) + rep( c(1,14,27,40), noYears+1)
      weeks <- subset(weeks, !is.na(match(weeks,quarterStarts)))
      weekIdx <- weeks - (firstweek-1)

      # get the right year for each week
      year <- weeks %/% 52 + startyear
      # function to define the quarter order
      quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
      # get the right number and order of quarter labels
      quarter <- sapply( (weeks-1) %/% 13 %% 4, quarterFunc)
    } else {   #If epochAsDate -- experimental functionality to handle ISO 8601
      date <- as.Date(x@epoch, origin="1970-01-01")
      years <- unique(as.numeric(formatDate(date,"%Y")))
      #Start of quarters in each year present in the data. 
      qStart <- as.Date(paste(rep(years,each=4), c("-01-01","-04-01","-07-01","-10-01"),sep=""))
      qName  <- rep(c("I","II","III","IV"), length.out=length(qStart))
      qIdx   <- qStart <= max(date)+10 & qStart >= min(date)-10
      qStart <- qStart[qIdx] ; qName <- qName[qIdx]
      #Find week in data closest to these dates
      weekIdx <- sapply(qStart, function(d) which.min(abs(as.numeric(date - d))))

      date <- date[weekIdx]
      #Year the ISO week belongs to
      year <- as.numeric(formatDate(date,"%G"))
      quarter <- qName
    }        
      
    #construct the computed axis labels -- add quarters if xaxis.units is requested
    if (xaxis.units) {
      labels.week <- paste(year,"\n\n",quarter,sep="")
    } else {
      labels.week <- paste(year,sep="")
    }

    axis( side=1,line=1,labels=FALSE,at=c(1,length(observed)),lwd.ticks=0)
    axis( at=weekIdx[which(quarter != "I")] , labels=labels.week[which(quarter != "I")] , side=1, line = 1 ,cex=cex)
    #Bigger tick marks at the first quarter
    at <- weekIdx[which(quarter == "I")]
    axis( at=at  , labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par()$tcl)
    #2nd axis
#    axis( side=2 ,cex=cex)
  } else { ##other frequency
    #A label at each unit
    myat.unit <- seq(firstweek,length.out=length(observed) )

    # get the right year order
    month <- (myat.unit-1) %% x@freq + 1
    year <- (myat.unit - 1) %/% x@freq + startyear
    #construct the computed axis labels -- add quarters if xaxis.units is requested
    if (xaxis.units) {
      mylabels.unit <- paste(year,"\n\n", (myat.unit-1) %% x@freq + 1,sep="")
    } else {
      mylabels.unit <- paste(year,sep="")
    }
    #Add axis
    axis( at=(1:length(observed))  , labels=NA, side=1, line = 1 ,cex=cex)
    axis( at=(1:length(observed))[month==1]  , labels=mylabels.unit[month==1] , side=1, line = 1 ,cex=cex)
#        axis( at=(1:length(observed)), labels=mylabels.unit, side=1, line = 1 ,cex=cex)
    #Bigger tick marks at the first unit
    at <- (1:length(observed))[(myat.unit - 1) %% x@freq == 0]
    axis( at=at  , labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par()$tcl)
    #2nd axis
#    axis( side=2 ,cex=cex)
  }
  invisible()
}

plot.sts.time.one <- function(x, k=1, domany=FALSE,ylim=NULL,xaxis.years=TRUE, axes=TRUE, xaxis.units=TRUE, epochsAsDate=x@epochAsDate, xlab="time", ylab="No. infected", main=NULL, type="s",lty=c(1,1,2),col=c(NA,1,4),lwd=c(1,1,1), outbreak.symbol = list(pch=3, col=3, cex=1),alarm.symbol=list(pch=24, col=2, cex=1),cex=1,legend.opts=list(x="top", legend=NULL,lty=NULL,pch=NULL,col=NULL),dx.upperbound=0.5,hookFunc=function() {},...) {

  #Extract slots -- depending on the algorithms: x@control$range
  observed   <- x@observed[,k]
  state      <- x@state[,k]
  alarm      <- x@alarm[,k]
  upperbound <- x@upperbound[,k]
  hasAlarm   <- all(!is.na(alarm))
  startyear <-  x@start[1]
  firstweek <-  x@start[2]
  method <-     x@control$name
  disease <-    x@control$data
  population <- x@populationFrac[,k]
  binaryTS <- x@multinomialTS

  if (binaryTS) {
    observed <- ifelse(population!=0,observed/population,0)
    upperbound <- ifelse(population!=0,upperbound/population,0)
    if (ylab == "No. infected") { ylab <- "Proportion infected" }
  }
  
   ##### Handle the NULL arguments ######################################
  if (is.null(main)) {
    #If no surveillance algorithm has been run
    if (length(x@control) != 0) {
      main = paste("Surveillance using ", as.character(method),sep="") 
    }
  }
  #No titles are drawn when more than one is plotted.
  if (domany) main = ""

  # control where the highest value is
  max <- max(c(observed,upperbound),na.rm=TRUE)
  
  #if ylim is not specified, give it a default value
  if(is.null(ylim) ){
    ylim <- c(-1/20*max, max)
  }

  # left/right help for constructing the columns
  dx.observed <- 0.5
  upperboundx <- (1:length(upperbound)) - (dx.observed - dx.upperbound)
  
  #Generate the matrices to plot (values,last value)
  xstuff <- cbind(c(upperboundx,length(observed) + min(1-(dx.observed - dx.upperbound),0.5)))
  ystuff <-cbind(c(upperbound,upperbound[length(observed) ]))

  #Plot the results 
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab=ylab,main=main,ylim=ylim,axes = !(xaxis.years),type=type,lty=lty[-c(1:2)],col=col[-c(1:2)],lwd=lwd[-c(1:2)],...)

  #This draws the polygons containing the number of counts (sep. by NA)
  i <- rep(1:length(observed),each=5)
  dx <- rep(dx.observed * c(-1,-1,1,1,NA), times=length(observed))
  x.points <- i + dx
  y.points <- as.vector(t(cbind(0, observed, observed, 0, NA)))
  polygon(x.points,y.points,col=col[1],border=col[2],lwd=lwd[1])

  #Draw upper bound once more in case the polygons are filled
  if (!is.na(col[1])) {
    lines(x=xstuff,y=ystuff,type=type,lty=lty[-c(1:2)],col=col[-c(1:2)],lwd=lwd[-c(1:2)],...)
  }
  
  #Draw outbreak symbols
  alarmIdx <- which(!is.na(alarm) & (alarm == 1))
  if (length(alarmIdx)>0) {
    matpoints( alarmIdx, rep(-1/40*ylim[2],length(alarmIdx)), pch=alarm.symbol$pch, col=alarm.symbol$col, cex= alarm.symbol$cex)
  }
  
  #Draw alarm symbols
  stateIdx <- which(state == 1)
  if (length(stateIdx)>0) {
    matpoints( stateIdx, rep(-1/20*ylim[2],length(stateIdx)), pch=outbreak.symbol$pch, col=outbreak.symbol$col,cex = outbreak.symbol$cex)
  }

  #Label x-axis 
  if(xaxis.years & axes) {
    addFormattedXAxis(x, epochsAsDate, observed, firstweek,xaxis.units,cex)
  }
  #Label y-axis
  if (axes) {
    axis( side=2 ,cex=cex)
  }

  if(!is.null(legend.opts)) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$x)) legend.opts$x <- "topleft"
    if (is.null(legend.opts$lty)) legend.opts$lty <- c(lty[1],lty[3],NA,NA)
    if (is.null(legend.opts$col)) legend.opts$col <- c(col[2],col[3],outbreak.symbol$col,alarm.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch <- c(NA,NA,outbreak.symbol$pch,alarm.symbol$pch)
    if (is.null(legend.opts$legend))
      legend.opts$legend <- c("Infected", "Threshold","Outbreak","Alarm" )
    #print(legend.opts)
    do.call("legend",legend.opts)
  }

  #Call hook function for user customized action
  environment(hookFunc) <- environment()
  hookFunc()

  invisible()
}


plot.sts.alarm <- function(x, lvl=rep(1,nrow(x)), ylim=NULL,xaxis.years=TRUE, xaxis.units=TRUE, epochsAsDate=x@epochAsDate, xlab="time", main=NULL, type="hhs",lty=c(1,1,2),col=c(1,1,4), outbreak.symbol = list(pch=3, col=3, cex=1),alarm.symbol=list(pch=24, col=2, cex=1),cex=1,cex.yaxis=1,...) {

  k <- 1
  #Extract slots -- depending on the algorithms: x@control$range
  observed   <- x@observed[,k]
  state      <- x@state[,k]
  alarm      <- x@alarm[,k]
  upperbound <- x@upperbound[,k]
  hasAlarm   <- all(!is.na(alarm))
  startyear <-  x@start[1]
  firstweek <-  x@start[2]
  method <-     x@control$name
  disease <-    x@control$data
  ylim <- c(0.5, ncol(x))
  
   ##### Handle the NULL arguments ######################################
  if (is.null(main)) {
    #If no surveillance algorithm has been run
    if (length(x@control) != 0) {
     # main = paste("Analysis of ", as.character(disease), " using ",
      main = paste("Surveillance using ", as.character(method),sep="") 
    }
  }
 
  # left/right help for constructing the columns
  dx.observed <- 0.5
  observedxl <- (1:length(observed))-dx.observed
  observedxr <- (1:length(observed))+dx.observed
  upperboundx <- (1:length(upperbound)) #-0.5
  
  # control where the highest value is
  max <- max(c(observed,upperbound),na.rm=TRUE)
        
  #if ylim is not specified
  if(is.null(ylim)){
    ylim <- c(-1/20*max, max)
  }


  #Generate the matrices to plot
  xstuff <- cbind(observedxl, observedxr, upperboundx)
  ystuff <-cbind(observed, observed, upperbound)
        

  #Plot the results using one Large plot call (we do this by modifying
  #the call). Move this into a special function!
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab="",main=main,ylim=ylim,axes = FALSE,type="n",lty=lty,col=col,...)

  #Label of x-axis 
  if(xaxis.years){
    addFormattedXAxis(x, epochsAsDate, observed, firstweek,xaxis.units,cex)
  }
  axis( side=2, at=1:ncol(x),cex.axis=cex.yaxis, labels=colnames(x),las=2)


  #Draw all alarms
  for (i in 1:nrow(x)) {
    idx <- (1:ncol(x))[x@alarm[i,] > 0]
    for (j in idx) {
      points(i,j,pch=3,col=lvl[j]+1)
    }
  }

  #Draw lines seperating the levels
  m <- c(-0.5,cumsum(as.numeric(table(lvl))))
  sapply(m, function(i) lines(c(0.5,nrow(x@alarm)+0.5),c(i+0.5,i+0.5),lwd=2))
  
  invisible()
}


#xaxis.years=TRUE,startyear = 2001, firstweek = 1, legend=TRUE
plot.sts.time <- function(x, type, method=x@control$name, disease=x@control$data,same.scale=TRUE,par.list=list(mfrow=magic.dim(nAreas),mar=par()$mar),...) {

  #Plot as one if type = time + unit 
  as.one=all(!is.na(pmatch(c("time","unit"),type[[3]] ))) & is.na(pmatch("|",type[[3]]))
  
  #Extract
  observed <- x@observed
  state <- x@state
  alarm <- x@alarm
  population <- x@populationFrac
  binaryTS <- x@multinomialTS
  
  #univariate timeseries ?
  if(is.vector(observed))
    observed <- matrix(observed,ncol=1)
  if(is.vector(state))
    state <- matrix(state,ncol=1)
  if(is.vector(alarm)) 
    alarm <- matrix(alarm,ncol=1)
  nAreas <- ncol(observed)

  if (binaryTS) {
    pi <-  ifelse(population == 0,  0,observed / population)
    if (identical(dim(x@upperbound), population)) {
      un <-  ifelse(population == 0, 0, x@upperbound / population)
    } else {
      un <- 0
    }
    max <-  max(max(pi),max(un),na.rm=TRUE)
  } else {
    max <-  max(max(observed),max(x@upperbound),na.rm=TRUE)
  }

  #Check empty arguments
  if (is.null(par.list[["mfrow",exact=TRUE]])) {
    par.list$mfrow <- magic.dim(nAreas)
  } 
  

  #multivariate time series
  if(nAreas > 1){
    #all areas in one plot
    if(as.one) {
      #### This is currently not supported
    } else {
      #set window size
      oldpar <- par()
      par(par.list)

      #All plots on same scale? If yes, then check if a scale
      #is already specified using the ylim argument
      args <- list(...)
      if(same.scale) {
        if (is.null(args$ylim)) {
          args$ylim <- c(-1/20*max, max)
        }
      } else {
        args$ylim <- NULL
      }
      
      #plot areas
      for (k in 1:nAreas) {
        #Changed call of plot.sts.time.one to invocation using "call"
        argsK <- merge.list(args,list("x"=x,"k"=k,"domany"=TRUE,"legend"=NULL))
        do.call("plot.sts.time.one",args=argsK)
        #Add title - do this using the cex.main size
        cex.main <- list(...)[["cex.main",exact=TRUE]]
        if (is.null(cex.main)) { cex.main <- 1 }
        mtext(colnames(observed)[k],line=-1.3,cex=cex.main)     
      }
      #reset graphical params
      #par(mfrow=c(1,1), mar=c(5, 4, 4, 2)+0.1)
      oldwarn <- options()$warn ; options(warn=-1) ; par(oldpar) ; options(warn=oldwarn)
    }
  } else {  #univariate time series
    plot.sts.time.one(x=x, domany=FALSE,...)
  }
  invisible()
}

##############################################################
# Animation
#
# opts.col - how to get colors??
# wait - waiting time (max count in a for loop) better way?
##############################################################

plot.sts.spacetime <- function(x,type,legend=NULL,opts.col=NULL,labels=TRUE,wait.ms=250,cex.lab=0.7,verbose=FALSE,dev.printer=NULL,...) {
  #Extract the mappoly
  if (nrow(x@map@data) == 0)
    stop("The sts object doesn't have a proper map entry.")
  map <- x@map
  maplim <- list(x=bbox(map)[1,],y=bbox(map)[2,])

  #Check colnames, otherwise no need to continue
  if (is.null(colnames(x@observed)))
    stop("The sts observed slot does not have any colnames to match with the shapefile.")

  #Check for color options
  if (is.null(opts.col)) {
    opts.col <- list(ncolors=100,use.color=TRUE)
  }
  #Check for legend options
  if (is.null(legend)) {
    legend <- list(dx=0.4,dy=0.04,x=maplim$x[1],y=maplim$y[1],once=TRUE)
  }

  #Process dev.printer options
  if (!is.null(dev.printer)) {
    #Device
    if (is.null(dev.printer$device)) dev.printer$device <- png
    #File extension
    if (is.null(dev.printer$extension)) dev.printer$extension <- ".png"
    #Width and height
    if (is.null(dev.printer$width)) dev.printer$width <- 640
    if (is.null(dev.printer$height)) dev.printer$height <- 480
  }
      

  #Extract the data
  o <- x@observed
  alarm <- x@alarm
  
  #Formula is of type "observed ~ 1|unit" (i.e. no time)
  aggregate <- type[[3]][[3]] == "unit"
  if (aggregate) {
    o <- t(as.matrix(apply(o,MARGIN=2,sum)))
    alarm <- t(as.matrix(apply(alarm,MARGIN=2,sum)))>0
  }
  
  #Number of time points
  maxt <- dim(o)[1]

  
  #Make colors using the hcl.colors function
  gyr <- hcl.colors(o,ncolors=length(o),use.color=TRUE)

  #Cut into specified number of colors
  o.cut <- matrix(as.numeric(cut(o,length(gyr$col))),nrow=nrow(o),ncol=ncol(o))
  o.col <- matrix(gyr$col[o.cut],ncol=ncol(o.cut))
  o.col[is.na(o.col)] <- gray(1)
  dimnames(o.col) <- dimnames(o)

  #Sort the o according to the names in the map
  region.id <- unlist(lapply(map@polygons,function(poly) poly@ID))
  o.col.id <- dimnames(o.col)[[2]]

  #Make the columns of o as in the map object
  o.col <- o.col[,pmatch(region.id,o.col.id),drop=FALSE]
  alarm.col <- alarm[,pmatch(region.id,o.col.id),drop=FALSE]

  #Screen processing
  screen.matrix <- matrix(c(0,1,0,1,0,1,0.8,1),2,4,byrow=TRUE)
  split.screen(screen.matrix)

  #Loop over all time slices
  for (t in 1:maxt) {
    #Status information
    if (verbose) {
      cat(paste("Processing slice",t,"of",maxt,"\n"))
    }
    
    #Clean screen (title area)
    screen(n=2)
    par(bg=gray(1))
    erase.screen()
    par(bg="transparent")

    #Plot the map on screen 1
    screen(n=1)
    plot(map,col=o.col[t,],xlab="",ylab="",...)
    #Indicate alarms as shaded overlays
    if (!all(is.na(alarm.col))) {
      #Plotting using density "NA" does not appear to work
      #anymore in the new sp versions
      alarm.col[is.na(alarm.col)] <- 0
      plot(map,dens=alarm.col*15,add=TRUE)
    }
    

    if (labels)
      #getSpPPolygonsLabptSlots is deprecated. Use coordinates method insteas
      text(coordinates(map), labels=as.character(region.id), cex.lab=cex.lab)
  
    if (!aggregate) { title(paste(t,"/",maxt,sep="")) }

    #In case a legend is requested
    if (!is.null(legend) && !(legend$once & t>1)  | (t==1)) {
      add.legend(legend,  maplim ,gyr)
    }

    #Is writing to files requested?
    if (!is.null(dev.printer)) {
      #Create filename
      fileName <- paste(dev.printer$name,"-",insert.zeroes(t,length=ceiling(log10(maxt))),
                        dev.printer$extension,sep="")
      cat("Creating ",fileName,"\n")
      #Save the current device using dev.print
     dev.print(dev.printer$device, file=fileName, width=dev.printer$width, height=dev.printer$height)
    }
    
    wait(wait.ms) 
  }
  close.screen(all.screens = TRUE)
}

###################################################################
#Helper function for waiting a specific amount of milliseconds
#implemented using proc.time. The waiting time is probably not 100%
#exact as R is an interpretated language.
#
# Params:
#  wait.ms - number of milliseconds to wait
###################################################################

wait <- function(wait.ms) {
  #Initialize
  start.time <- proc.time()[3]*1000
  ellapsed <- proc.time()[3]*1000 - start.time

  #Loop as long as required.
  while (ellapsed < wait.ms) {
    ellapsed <- proc.time()[3]*1000 - start.time
  }
}

#Helper function -- use colors from the vcd package
hcl.colors <- function(x, ncolors=100, use.color=TRUE)
{
    GYR <- if (require("colorspace")) { # the Zeil-ice colors 
        colorspace::heat_hcl(ncolors, h=c(0,120),
                             c=if (use.color) c(90,30) else c(0,0),
                             l=c(50,90), power=c(0.75, 1.2))
    } else {
        if (use.color) heat.colors(ncolors) else grey.colors(ncolors)
    }
    list(col=rev(GYR), min=0, max=max(x), trans=base::identity)
}

#Helper function for plot.sts.spacetime
add.legend <- function(legend,maplim,theColors) {
  #Preproc
  dy <- diff(maplim$y) * legend$dy
  dx <- diff(maplim$x) * legend$dx
    
  #Add legend -- i.e. a slider
  xlu <- xlo <- legend$x
  xru <- xro <- xlu + dx 
  yru <- ylu <- legend$y
  yro <- ylo <- yru + dy 

  
  step <- (xru - xlu)/length(theColors$col)
  for (i in 0:(length(theColors$col) - 1)) {
    polygon(c(xlo + step * i, xlo + step * (i + 1), 
              xlu + step * (i + 1), xlu + step * i), c(ylo, 
                                                       yro, yru, ylu), col = theColors$col[i + 1], 
            border = theColors$col[i + 1])
  }
  
  
  #Write info about min and max on the slider.
  black <- grey(0)
  lines(c(xlo, xro, xru, xlu, xlo), c(ylo, yro, yru, ylu, ylo), col =   black)

  #Transformation function for data values, either
  #exp or identical
  trans <- theColors$trans

  text(xlu, ylu - 0.5*dy, formatC(trans(theColors$min)), cex = 1, col = black,adj=c(0,1))
  text(xru, yru - 0.5*dy, formatC(trans(theColors$max)), cex = 1, col = black,adj=c(1,1))
}

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

  ## FIXME: if map slot is not NULL, check that all colnames(object@observed)
  ## are in row.names(object@map) (ideally, these are identical sets)
  
  if(is.null( retval)) return ( TRUE )
  else return ( retval )
})



######################################################################
# Insert leading zeros so integers obtain a fixed length. Good when
# making filenames so they are all of same length ideal for sorting.
#
# Parameters
#  x - An integer
######################################################################

insert.zeroes<- function(x,length=3) {
  for (i in 1:(length-1)) {
    if (x<10^i) return(paste(paste(rep(0,length-i),collapse=""),x,sep=""))
  }
  #If x has more digits than length then just return x
  return(paste(x))
}


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

  cat("\nmap:\n")
  print(object@map@data$SNAME)

  cat("\nhead of neighbourhood:\n")
  print( head(object@neighbourhood,n))
} )

######################################################################
#Method to convert sts object to a data frame suitable for regression
#Params:
# row.names - from generic R function
# optional  - from generic R function
# freqByWeek -- if TRUE use information in week (supposed to be Dates)
#               to freq (e.g. used for regression model)
######################################################################

setMethod("as.data.frame", signature(x="sts"), function(x,row.names = NULL, optional = FALSE, ...) {
  #Convert object to data frame and give names
  res <- data.frame("observed"=x@observed, "epoch"=x@epoch, "state"=x@state, "alarm"=x@alarm,"population"=x@populationFrac)

  if (ncol(x) > 1) {
    colnames(res) <-  c(paste("observed.",colnames(x@observed),sep=""),"epoch",
                        paste("state.",colnames(x@observed),sep=""),
                        paste("alarm.",colnames(x@observed),sep=""),
                        paste("population.",colnames(x@observed),sep=""))
  } else {
      colnames(res) <-  c("observed","epoch","state","alarm","population")
  }
  
  #Add a column denoting the number of week
  if (x@epochAsDate) {
    #Convert to date
    date <- as.Date(x@epoch, origin="1970-01-01")
    epochStr <- switch( as.character(x@freq), 
                       "12" = "%m",
                       "52" =  "%V",
                       "365" = "%j")
                       
    #Find out how many epochs there are each year
    years <- unique(as.numeric(formatDate(date,"%Y")))
    dummyDates <- as.Date(paste(rep(years,each=6),"-12-",26:31,sep=""))
    maxEpoch <- tapply( as.numeric(formatDate(dummyDates, epochStr)), rep(years,each=6), max)
    #Assign this to result
    res$freq <- maxEpoch[pmatch(formatDate(date,"%Y"),names(maxEpoch),duplicates.ok=TRUE)]
    res$epochInPeriod <- as.numeric(formatDate(date,epochStr)) / res$freq
  } else {
    #Otherwise just replicate the fixed frequency
    res$freq <- x@freq
    res$epochInPeriod <- x@epoch %% res$freq
  }
  
  return(res)
})

