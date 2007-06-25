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
init.sts <- function(.Object, week, start=c(2000,1), freq=52, observed, state, map=NULL, neighbourhood=NULL, populationFrac=NULL,alarm=NULL,upperbound=NULL, control=NULL) {
  #Name handling
  namesObs <-colnames(observed)
  namesState <- colnames(observed)

  #check number of columns of observed and state
  nAreas <- ncol(observed)
  nObs <- nrow(observed)
  if(ncol(observed) != ncol(state)){
    #if there is only one state-vector for more than one area, repeat it
    if(ncol(state)==1)
      state <- ts(matrix(rep(state,nAreas),ncol=nAreas,byrow=FALSE),freq=frequency(observed))
    else{ 
      cat('wrong dimensions of observed and state \n')
      return(NULL)
    }
  }
  
  #check neighbourhood matrix
  if(!is.null(neighbourhood) & (any(dim(neighbourhood) != nAreas))) {
    cat('wrong dimensions of neighbourhood matrix \n')
    return(NULL)
  }

  #popFrac
  if (is.null(populationFrac)) {
    populationFrac <- matrix(1/nAreas,nrow=nObs,ncol=nAreas)
  }
  if (nAreas ==1){
    populationFrac <- matrix(1,nrow=nObs, ncol=1)
  }
  
  #labels for observed and state
  if(is.null(namesObs)){
    namesObs <- paste("observed", 1:nAreas, sep="")       
    namesState <- paste("state", 1:nAreas, sep="")  
  }
 
  dimnames(observed) <- list(NULL,namesObs)
  dimnames(state) <- list(NULL,namesState)

  if (is.null(neighbourhood))
    neighbourhood <- matrix(NA,nrow=ncol(observed),ncol=ncol(observed))
  if (is.null(alarm)) 
    alarm      <- matrix(NA,nrow=dim(observed)[1],ncol=dim(observed)[2])
  if (is.null(upperbound))
    upperbound <- matrix(NA,nrow=dim(observed)[1],ncol=dim(observed)[2])

  ##Assign everything
  .Object@week <- week

  if (length(start) == 2) {
    .Object@start <- start
  } else {
    stop("start must be a vector of length two denoting (year, week/month/idx)")
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
  
  return(.Object)
}

###########################################################################
# Initialization -- two modes possible: full or just disProg, freq and map
###########################################################################

#Full
setMethod("initialize", "sts", init.sts)

#Partial -- use a disProg object as start and convert it.
disProg2sts <- function(disProgObj, map=NULL) {
  sts <- new("sts", week=disProgObj$week, start=disProgObj$start, freq=disProgObj$freq, observed=disProgObj$observed, state = disProgObj$state, map=map, neighbourhood=disProgObj$neighbourhood, populationFrac=disProgObj$populationFrac,alarm=disProgObj$alarm,upperbound=disProgObj$upperbound)
  return(sts)
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
    x@week <- 1:m
    
    x@observed <- as.matrix(aggregate(x@observed,by=list(new),sum)[,-1])
    x@state <- as.matrix(aggregate(x@state,by=list(new),sum)[,-1])>0
    x@alarm <- as.matrix(aggregate(x@alarm,by=list(new),sum)[,-1])
    x@upperbound <- as.matrix(aggregate(x@upperbound,by=list(new),sum)[,-1])
    x@populationFrac <- as.matrix(aggregate(x@populationFrac,by=list(new),sum)[,-1])
  }
  if (by == "unit") {
    #Aggregate units
    x@observed <- as.matrix(apply(x@observed, MARGIN=1, sum))
    x@state <- as.matrix(apply(x@state, MARGIN=1, sum))>0
    x@alarm <- as.matrix(apply(x@alarm, MARGIN=1, sum))>0
    #There is no clever way to aggregate the upperbounds
    x@upperbound <- matrix(NA,ncol=ncol(x@alarm),nrow=nrow(x@alarm))
    x@populationFrac <- as.matrix(apply(x@populationFrac, MARGIN=1, sum))>0
    x@neighbourhood <- matrix(1,nrow=1,ncol=1)
  }

  #validObject(x) #just a check

  return(x)
})


#####################################################################
# Misc access functions
####################################################################

setMethod("nrow", "sts", function(x) return(nrow(x@observed)))
setMethod("ncol", "sts", function(x) return(ncol(x@observed)))
setMethod("dim", "sts", function(x) return(dim(x@observed)))
setMethod("colnames", signature=c(x="sts",do.NULL="missing",prefix="missing"), function(x,do.NULL, prefix) return(colnames(x@observed)))

#####################################################################
#[-method for accessing the observed, alarm, etc. objects
#####################################################################

setMethod("[", "sts", function(x, i, j, ..., drop) {
  #default value for i and j
  if(missing(i)) {i <- min(1,nrow(x@observed)):nrow(x@observed)}
  if(missing(j)) {j <- min(1,ncol(x@observed)):ncol(x@observed)}

  x@week <- x@week[i]
  x@observed <- x@observed[i,j,drop=FALSE]
  x@state <- x@state[i,j,drop=FALSE]
  x@alarm <- x@alarm[i,j,drop=FALSE]
  x@populationFrac <- x@populationFrac[i,j,drop=FALSE]
  x@upperbound <- x@upperbound[i,j,drop=FALSE]

  #Neighbourhood matrix
  x@neighbourhood <- x@neighbourhood[j,j,drop=FALSE]
  
  #Fix the corresponding start entry
  start <- x@start
  new.sampleNo <- start[2] + min(i) - 1
  start.year <- start[1] + (new.sampleNo - 1) %/% x@freq 
  start.sampleNo <- (new.sampleNo - 1) %% x@freq + 1
  x@start <- c(start.year,start.sampleNo)

  #Save time by not allocating a new object
  #res <- new("sts",week=week, freq=x@freq, start=start,observed=observed,state=state,alarm=alarm,upperbound=upperbound,neighbourhood=neighbourhood,populationFrac=populationFrac,map=x@map,control=x@control)
    
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
  obsOk <- type[[2]] == "observed"
  map   <- (length(type[[3]])==3) && (type[[3]][[1]] == "|") && (type[[3]][[2]] == "1")
  time  <- pmatch("time",type[[3]]) > 0

  #Valid formula?
  valid <- lapply(as.list(type[[3]]),function(i) is.na(pmatch(i,c("1","unit","|","time","*","+"))))
  valid <- all(!unlist(valid))

  #No unit dimenstion?
  justTime <- type[[3]] == "time"
  
  if (!obsOk | !valid) {
    stop("Not a valid plot type.")
  }


  #space-time plots
  if (map) {
    plot.sts.spacetime(x,type,...)
    return(invisible())
  }
  #time plots
  if (time) {
    #In case observed ~ time, the units are aggregated
    plot.sts.time( if(justTime) aggregate(x,by="unit") else x,type,...)
    return(invisible())
  }
})


##########################################################################
# Plot functions
##########################################################################

plot.sts.time.one <- function(x, k=1, domany=FALSE,ylim=NULL,xaxis.years=TRUE, xaxis.units=TRUE, xlab="time", ylab="No. infected", main=NULL, type="hhs",lty=c(1,1,2),col=c(1,1,4), outbreak.symbol = list(pch=3, col=3, cex=1),alarm.symbol=list(pch=24, col=2, cex=1),cex=1,legend.opts=list(x="top", legend=NULL,lty=NULL,pch=NULL,col=NULL),...) {

  #Extract slots -- depending on the algorithms: x@control$range
  observed   <- x@observed[,k]
  state      <- x@state[,k]
  alarm      <- x@alarm[,k]
  upperbound <- x@upperbound[,k]
  hasAlarm   <- all(!is.na(alarm))
  startyear <- x@start[1]
  firstweek <- x@start[2]
  method <- x@control$name
  disease <- x@control$data
  
   ##### Handle the NULL arguments ######################################
  if (is.null(main)) {
    #If no surveillance algorithm has been run
    if (length(x@control) != 0) {
     # main = paste("Analysis of ", as.character(disease), " using ",
      main = paste("Surveillance using ", as.character(method),sep="") 
    }
  }
  #No titles are drawn when more than one is plotted.
  if (domany) main = ""

 
  # width of the column
  tab <- 0.5

  # left/right help for constructing the columns
  observedxl <- (1:length(observed))-tab
  observedxr <- (1:length(observed))+tab
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
  #the call). 
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab=ylab,main=main,ylim=ylim,axes = !(xaxis.years),type=type,lty=lty,col=col,...)

  #Aggregated data types -- this is not supported atm
  #if (!is.null(survResObj$aggr)) {
  #  points(upperboundx+tab,survResObj$aggr,col=1)
  #}
  
  for(i in 1:length(observed)){
    matlines( c(i-tab, i+tab), c(observed[i],observed[i]), col=col[1] )
    if (!is.na(alarm[i]) && (alarm[i] == 1))
      matpoints( i, -1/40*max, pch=alarm.symbol$pch, col=alarm.symbol$col, cex= alarm.symbol$cex)
    if (state[i] == 1)
      matpoints( i, -1/20*max, pch=outbreak.symbol$pch, col=outbreak.symbol$col,cex = outbreak.symbol$cex)
  }

  # check where to place the legend. If the left upper side is free place it there
  if (max * 2/3 >= max(
                max(observed[1:floor(1/4 * length(observed))]),
                max(upperbound[1:floor(1/4 * length(upperbound))]),na.rm=TRUE
                )) {
    xlegpos <- 0
  }
  #Label of x-axis 
  if(xaxis.years){
    if (x@freq ==52) {
      # get the number of quarters lying in range for getting the year and quarter order
      myat.week <- seq(ceiling((52-firstweek+1)/13) * 13 + 1, length(observed)+(floor((52-firstweek + 1)/13) * 13 +1), by=13)
      # get the right year order
      year <- (myat.week - 52) %/% 52 + startyear
      # function to define the quarter order
      quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
      # get the right number and order of quarter labels
      quarter <- sapply( (myat.week-1) %/% 13 %% 4, quarterFunc)
      # get the positions for the axis labels
      myat.week <- myat.week - (52 - firstweek + 1)

      #construct the computed axis labels -- add quarters if xaxis.units is requested
      if (xaxis.units) {
        mylabels.week <- paste(year,"\n\n",quarter,sep="")
      } else {
        mylabels.week <- paste(year,sep="")
      }
    
      axis( at=myat.week , labels=mylabels.week , side=1, line = 1 ,cex=cex)
      axis( side=2 ,cex=cex)
    } else { ##other frequency
      #A label at each unit
      myat.unit <- seq(firstweek,length.out=length(observed) )

      # get the right year order
      year <- (myat.unit - 1) %/% x@freq + startyear
      #construct the computed axis labels -- add quarters if xaxis.units is requested
      if (xaxis.units) {
        mylabels.unit <- paste(year,"\n\n", (myat.unit-1) %% x@freq + 1,sep="")
      } else {
        mylabels.unit <- paste(year,sep="")
      }
      #Add axis
      axis( at=1:length(observed)  , labels=mylabels.unit , side=1, line = 1 ,cex=cex)
      #Bigger tick marks at the first unit
      at <- (1:length(observed))[(myat.unit - 1) %% x@freq == 0]
      axis( at=at  , labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par()$tcl)
      #2nd axis
      axis( side=2 ,cex=cex)
    }
  }

  if(!is.null(legend.opts)) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$lty)) legend.opts$lty = c(lty[1],lty[3],NA,NA)
    if (is.null(legend.opts$col)) legend.opts$col = c(col[1],col[3],outbreak.symbol$col,alarm.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch = c(NA,NA,outbreak.symbol$pch,alarm.symbol$pch)
    if (is.null(legend.opts$legend))
      legend.opts$legend = c("Infected", "Threshold","Outbreak","Alarm" )
    
    do.call("legend",legend.opts)
  }

  invisible()
}


#xaxis.years=TRUE,startyear = 2001, firstweek = 1, legend=TRUE
plot.sts.time <- function(x, type, method=x@control$name, disease=x@control$data,same.scale=TRUE,...) {

  #Plot as one if type = time + unit 
  as.one=all(!is.na(pmatch(c("time","unit"),type[[3]] ))) & is.na(pmatch("|",type[[3]]))
  
  #Extract
  observed <- x@observed
  state <- x@state
  alarm <- x@alarm

  #univariate timeseries ?
  if(is.vector(observed))
    observed <- matrix(observed,ncol=1)
  if(is.vector(state))
    state <- matrix(state,ncol=1)
  if(is.vector(alarm)) 
    alarm <- matrix(alarm,ncol=1)
  nAreas <- ncol(observed)
  
  max <-  max(max(observed),max(x@upperbound),na.rm=TRUE)

  #multivariate time series
  if(nAreas > 1){
    #all areas in one plot
    if(as.one) {
      #### This is currently not supported
    } else {
      #set window size     
      par(mfrow=magic.dim(nAreas),mar=c(2,1,1,1))
      
      if(same.scale) {
        ylim <- c(-1/20*max, max)
      } else {
        ylim <- NULL
      }
      
      #plot areas
      k <- 1:nAreas
      sapply(k, function(k) {
        plot.sts.time.one(x, k=k, domany=TRUE, ylim=ylim, legend=NULL, ... )   
        mtext(colnames(observed)[k],line=-1.3)     
      })
      #reset graphical params
      par(mfrow=c(1,1), mar=c(5, 4, 4, 2)+0.1)
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

  #Sort the o xected according to the names in the map xect
  region.id <- unlist(lapply(map@polygons,function(poly) poly@ID))
  o.col.id <- dimnames(o.col)[[2]]

  #Make the columns of o as in the map object
  o.col <- o.col[,pmatch(region.id,o.col.id),drop=FALSE]

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
    screen(2)
    par(bg=gray(1))
    erase.screen()
    par(bg="transparent")

    #Plot the map on screen 1
    screen(1)
    plot(map,col=o.col[t,],xlab="",ylab="",...)
    plot(map,dens=alarm*15,add=TRUE)
    

    if (labels)
      text(getSpPPolygonsLabptSlots(map), labels=as.character(region.id), cex.lab=cex.lab)
  
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
hcl.colors <- function(x,ncolors=100,use.color=TRUE) {
  if (use.color) {
    #The Zeil-ice colors 
    GYR <- rev(heat_hcl(ncolors, h=c(0,120), c=c(90,30), l=c(50,90), power=c(0.75, 1.2)))
  } else {
    #Sanity check
    GYR <- rev(heat_hcl(ncolors, h=c(0,120), c=0, l=c(50,90), power=c(0.75, 1.2)))
  }
  return(list(col=GYR,min=0,max=max(x), trans=function(x) return(x)))
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
  #cat( "length(week):\t", length(object@week),"\n" )
  cat( "freq:\t\t", object@freq,"\n" )
  cat( "start:\t\t",object@start,"\n" )
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



