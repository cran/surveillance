################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Time series plot for sts-objects
###
### Copyright (C) 2007-2013 Michael Hoehle
### $Revision: 918 $
### $Date: 2014-05-10 23:16:41 +0200 (Sat, 10 May 2014) $
################################################################################


######################################################################
# stsplot_time sets the scene and calls stsplot_time1 for each unit
######################################################################

stsplot_time <- function(x, method=x@control$name, disease=x@control$data,
                         as.one=FALSE, same.scale=TRUE, par.list=list(), ...)
{
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

  #If no "mar" or "mfrow" argument in par.list add default values.
  if (is.null(par.list[["mar",exact=TRUE]])) { par.list$mar <- c(5,4,1,1)}
  if (is.null(par.list[["mfrow",exact=TRUE]])) { par.list$mfrow=magic.dim(nAreas)}
  
  #multivariate time series
  if(nAreas > 1){
    if(as.one) { # all areas in one plot
      stop("this type of plot is currently not implemented")
    } else {
      #set window size
      oldpar <- par(par.list)

      #All plots on same scale? If yes, then check if a scale
      #is already specified using the ylim argument
      args <- list(...)
      if(same.scale) {
        if (is.null(args$ylim)) {
            max <- if (binaryTS) {
                max(ifelse(population == 0, 0,
                           pmax(observed,x@upperbound,na.rm=TRUE)/population))
            } else max(observed,x@upperbound,na.rm=TRUE)
            args$ylim <- c(-1/20*max, max)
        }
      } else {
        args$ylim <- NULL
      }
      
      #plot areas
      for (k in 1:nAreas) {
        argsK <- modifyList(args, list(x=x, k=k, main="", legend.opts=NULL),
                            keep.null = TRUE)
        do.call("stsplot_time1",args=argsK)
        title(main=colnames(observed)[k], line=-1)
      }
      
      #reset graphical params
      par(oldpar)
    }
  } else {  #univariate time series
    par.list$mfrow <- NULL #no mf formatting..
    oldpar <- par(par.list)
    stsplot_time1(x=x, ...)
    par(oldpar)
  }
  invisible()
}


### work-horse which produces the time series plot with formatted x-axis

stsplot_time1 <- function(
    x, k=1, ylim=NULL, axes=TRUE, xaxis.tickFreq=list("%Q"=atChange),
    xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="%G\n\n%OQ",
    epochsAsDate=x@epochAsDate, xlab="time", ylab="No. infected", main=NULL,
    type="s", lty=c(1,1,2), col=c(NA,1,4), lwd=c(1,1,1),
    outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1),
    alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
    legend.opts=list(x="top", legend=NULL, lty=NULL, pch=NULL, col=NULL),
    dx.upperbound=0L, hookFunc=function(){}, ...)
{

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

  #Control what axis style is used
  xaxis.dates <- !is.null(xaxis.labelFormat) 

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
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab=ylab,main=main,ylim=ylim,axes = !(xaxis.dates),type=type,lty=lty[-c(1:2)],col=col[-c(1:2)],lwd=lwd[-c(1:2)],...)

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
    matpoints( alarmIdx, rep(-1/40*ylim[2],length(alarmIdx)), pch=alarm.symbol$pch, col=alarm.symbol$col, cex= alarm.symbol$cex, lwd=alarm.symbol$lwd)
  }
  
  #Draw alarm symbols
  stateIdx <- which(state == 1)
  if (length(stateIdx)>0) {
    matpoints( stateIdx, rep(-1/20*ylim[2],length(stateIdx)), pch=outbreak.symbol$pch, col=outbreak.symbol$col,cex = outbreak.symbol$cex,lwd=outbreak.symbol$lwd)
  }

  #Label x-axis 
  if(xaxis.dates & axes) {
    addFormattedXAxis(x, epochsAsDate, observed, firstweek,xaxis.labelFormat,xaxis.tickFreq,xaxis.labelFreq,...)
  }
  #Label y-axis
  if (axes) {
    axis( side=2 ,...)#cex=cex, cex.axis=cex.axis)
  }

  if(is.list(legend.opts)) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$x)) legend.opts$x <- "topleft"
    if (is.null(legend.opts$lty)) legend.opts$lty <- c(lty[1],lty[3],NA,NA)
    if (is.null(legend.opts$col)) legend.opts$col <- c(col[2],col[3],outbreak.symbol$col,alarm.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch <- c(NA,NA,outbreak.symbol$pch,alarm.symbol$pch)
    if (is.null(legend.opts$legend)) {
      legend.opts$legend <- c("Infected", "Threshold","Outbreak","Alarm" )
    }
    #Make the legend
    do.call("legend",legend.opts)
  }

  #Call hook function for user customized action
  environment(hookFunc) <- environment()
  hookFunc()

  invisible()
}


##############
### alarm plot
##############

stsplot_alarm <- function(
    x, lvl=rep(1,nrow(x)), ylim=NULL,
    xaxis.tickFreq=list("%Q"=atChange),
    xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="%G\n\n%OQ",
    epochsAsDate=x@epochAsDate, xlab="time", main=NULL,
    type="hhs", lty=c(1,1,2), col=c(1,1,4),
    outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1),
    alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
    cex=1, cex.yaxis=1, ...)
{

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
 
  #Control what axis style is used
  xaxis.dates <- !is.null(xaxis.labelFormat) 

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
  if(xaxis.dates){
    addFormattedXAxis(x, epochsAsDate, observed, firstweek, xaxis.labelFormat,xaxis.tickFreq,xaxis.labelFreq,...)
  }
  axis( side=2, at=1:ncol(x),cex.axis=cex.yaxis, labels=colnames(x),las=2)


  #Draw all alarms
  for (i in 1:nrow(x)) {
    idx <- (1:ncol(x))[x@alarm[i,] > 0]
    for (j in idx) {
      points(i,j,pch=alarm.symbol$pch,col=alarm.symbol$col[lvl[j]],cex=alarm.symbol$cex,lwd=alarm.symbol$lwd)
    }
  }

  #Draw lines seperating the levels
  m <- c(-0.5,cumsum(as.numeric(table(lvl))))
  sapply(m, function(i) lines(c(0.5,nrow(x@alarm)+0.5),c(i+0.5,i+0.5),lwd=2))
  
  invisible()
}



#####################################
### Utilities to set up the time axis
#####################################

#Every unit change
atChange <- function(x,xm1) {
    which(diff(c(xm1,x)) != 0)
}
#Median index of factor
atMedian <- function(x,xm1) {
    as.integer(tapply(seq_along(x), INDEX=x, quantile, prob=0.5,type=3))
}
#Only every second unit change
at2ndChange <- function(x,xm1) {
    idxAtChange <- atChange(x,xm1)
    idxAtChange[seq(idxAtChange) %% 2 == 1]
}

#Helper function to format the x-axis of the time series
addFormattedXAxis <- function(x, epochsAsDate, observed, firstweek, xaxis.labelFormat, xaxis.tickFreq,xaxis.labelFreq,...) {
  
  #Old style if there are no Date objects
  if (!epochsAsDate) {
    #Declare commonly used variables.
    startyear <-  x@start[1]
    
    if (x@freq ==52) {  #Weekly epochs are the most supported 
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
      quarterFunc <- function(i) { switch(i+1,"I","II","III","IV") } #nicer:as.roman, but changes class.
      # get the right number and order of quarter labels
      quarter <- sapply( (weeks-1) %/% 13 %% 4, quarterFunc)
      #Computed axis labels -- add quarters (this is the old style)
      labels.week <- paste(year,"\n\n",quarter,sep="")
      #Make the line. Use lwd.ticks to get full line but no marks.
      axis( side=1,labels=FALSE,at=c(1,length(observed)),lwd.ticks=0,line=1,...) 
      axis( at=weekIdx[which(quarter != "I")] , labels=labels.week[which(quarter != "I")] , side=1, line = 1 ,...)       
      #Bigger tick marks at the first quarter (i.e. change of the year)
      at <- weekIdx[which(quarter == "I")]
      axis( at=at, labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par()$tcl)
      
    } else { ##other frequency (not really supported)
      #A label at each unit
      myat.unit <- seq(firstweek,length.out=length(observed) )
      
      # get the right year order
      month <- (myat.unit-1) %% x@freq + 1
      year <- (myat.unit - 1) %/% x@freq + startyear
      #construct the computed axis labels -- add quarters if xaxis.units is requested
      mylabels.unit <- paste(year,"\n\n", (myat.unit-1) %% x@freq + 1,sep="")
      
      #Add axis
      axis( at=(1:length(observed)), labels=NA, side=1, line = 1, ...) 
      axis( at=(1:length(observed))[month==1], labels=mylabels.unit[month==1] , side=1, line = 1 ,...)
      #Bigger tick marks at the first unit
      at <- (1:length(observed))[(myat.unit - 1) %% x@freq == 0]
      axis( at=at, labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par()$tcl)
    } 
  } else {   
    ################################################################
    #epochAsDate -- experimental functionality to handle ISO 8601
    ################################################################
    
    dates <- epoch(x)
    #make one which has one extra element at beginning with same spacing
    datesOneBefore <- c(dates[1]-(dates[2]-dates[1]),dates)
    
    #Make the line. Use lwd.ticks to get full line but no marks.
    axis( side=1,labels=FALSE,at=c(1,length(observed)),lwd.ticks=0,...) 
    
    ###Make the ticks (depending on the selected level).###
    tcl <- par()$tcl
    tickFactors <- surveillance.options("stsTickFactors")
    
    #Loop over all pairs in the xaxis.tickFreq list
     for (i in seq(xaxis.tickFreq)) {
      format <- names(xaxis.tickFreq)[i]
      xm1x <- as.numeric(formatDate(datesOneBefore,format))                         
      idx <- xaxis.tickFreq[[i]](x=xm1x[-1],xm1=xm1x[1])
      
      #Find tick size by table lookup
      tclFactor <- tickFactors[pmatch(format, names(tickFactors))]
      if (is.na(tclFactor)) { 
        warning("No tcl factor found for",format,". Setting it at 1") 
        tclFactor <- 1
      }
      axis(1,at=idx, labels=NA,tcl=tclFactor*tcl,...)
    }  
    
    ###Make the labels (depending on the selected level)###
    if (!is.null(xaxis.labelFormat)) {
      labelIdx <- NULL
      for (i in seq(xaxis.labelFreq)) {
        format <- names(xaxis.labelFreq)[i]
        xm1x <- as.numeric(formatDate(datesOneBefore,format))                         
        labelIdx <- c(labelIdx,xaxis.labelFreq[[i]](x=xm1x[-1],xm1=xm1x[1]))
      }
      
      #Format labels (if any) for the requested subset
      if (length(labelIdx)>0) {
          labels <- rep(NA,nrow(x))
          labels[labelIdx] <- formatDate(epoch(x)[labelIdx],xaxis.labelFormat)
          axis(1,at=1:nrow(x), labels=labels,tick=FALSE,...)
      }
    }
  }#end epochAsDate
  #Done
  invisible()
}
