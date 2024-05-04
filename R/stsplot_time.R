################################################################################
### Time series plot for sts-objects
###
### Copyright (C) 2007-2014 Michael Hoehle, 2013-2016,2021-2022 Sebastian Meyer
###
### This file is part of the R package "surveillance",
### free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at https://www.R-project.org/Licenses/.
################################################################################


######################################################################
# stsplot_time() sets the scene and calls either stsplot_time_as1()
# or stsplot_time1() for each unit
######################################################################

stsplot_time <- function(x, units = NULL,
                         as.one = FALSE, same.scale = TRUE,
                         par.list = list(), ...)
{
  observed <- x@observed
  if (is.null(units)) # plot all units
    units <- seq_len(ncol(observed))
  nUnits <- length(units)

  #graphical parameters
  if (is.list(par.list)) {
    if (nUnits > 1 && !as.one) {
      par.list <- modifyList( #default: reduced margins and mfrow panels
        list(mar = c(5,4,1,1), mfrow = magic.dim(nUnits)),
        par.list)
    } else {
      par.list$mfrow <- NULL #no mf formatting..
    }
    if (length(par.list) > 0) {
      oldpar <- par(par.list)
      on.exit(par(oldpar))
    }
  }

  if (nUnits == 1L) { # a single time-series plot
    stsplot_time1(x = x, k = units, ...)
  } else { # multiple time series
    if (as.one) { # all time series in one plot
      stsplot_time_as1(x, units = units, ...)
    } else { # each time series in a separate plot
      args <- list(...)
      if(same.scale) { # compute suitable ylim if not specified
        if (is.null(args[["ylim"]])) {
            ymax <- if (x@multinomialTS) {
                max(0, pmax(observed,x@upperbound,na.rm=TRUE)/x@populationFrac,
                    na.rm=TRUE)
            } else {
                max(observed,x@upperbound,na.rm=TRUE)
            }
            args$ylim <- c(-1/20*ymax, ymax)
        }
      } else {
        args$ylim <- NULL
      }

      #plot areas
      for (k in units) {
        argsK <- modifyList(args, list(x=x, k=k, main="", legend.opts=NULL),
                            keep.null = TRUE)
        do.call("stsplot_time1",args=argsK)
        title(main=if (is.character(k)) k else colnames(observed)[k], line=-1)
      }
    }
  }

  invisible()
}


## a simple matplot of observed counts from all/selected units, with a legend

stsplot_time_as1 <- function (x, units = NULL,
    type = "l", lty = 1:5, lwd = 1, col = 1:6,
    epochsAsDate = x@epochAsDate, xaxis.tickFreq = list("%Q"=atChange),
    xaxis.labelFreq = xaxis.tickFreq, xaxis.labelFormat = "%G\n\n%OQ",
    xlab = "time", ylab = "No. infected", legend.opts = list(), ...)
{
    observed <- x@observed
    if (x@multinomialTS) {
        observed <- ifelse(x@populationFrac != 0, observed/x@populationFrac, 0)
    }
    if (!is.null(units))
        observed <- observed[, units, drop = FALSE]

    ## basic plot
    opar <- par(bty = "n", xaxt = "n")  # a formatted time axis is added below
    matplot(observed, type = type, lty = lty, lwd = lwd, col = col,
            xlab = xlab, ylab = ylab, ...)
    par(opar)

    ## add time axis
    xaxis.line <- !epochsAsDate || grepl("\n", xaxis.labelFormat)
    addFormattedXAxis(x = x, epochsAsDate = epochsAsDate,
        xaxis.tickFreq = xaxis.tickFreq, xaxis.labelFreq = xaxis.labelFreq,
        xaxis.labelFormat = xaxis.labelFormat)  # line = 1

    ## add legend
    if (is.list(legend.opts)) {
        legend.opts <- modifyList(
            list(x = "top", legend = colnames(observed), lty = lty, lwd = lwd,
                 col = col, ncol = magic.dim(ncol(observed))[2L],
                 bty = "n"),
            legend.opts)
        do.call("legend", legend.opts)
    }

    invisible()
}


### work-horse which produces a single time series plot with formatted x-axis

stsplot_time1 <- function(
    x, k=1, ylim=NULL, axes=TRUE, xaxis.tickFreq=list("%Q"=atChange),
    xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="%G\n\n%OQ",
    epochsAsDate=x@epochAsDate, xlab="time", ylab="No. infected", main=NULL,
    type="s", lty=c(1,1,2), col=c(NA,1,4), lwd=c(1,1,1),
    outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1),
    alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
    legend.opts=list(),
    dx.upperbound=0L, hookFunc=function(){}, .hookFuncInheritance=function() {}, ...)
{
  stopifnot(length(k) == 1, is.character(k) || k != 0)

  #Extract slots -- depending on the algorithms: x@control$range
  observed   <- x@observed[,k]
  state      <- x@state[,k]
  alarm      <- x@alarm[,k]
  upperbound <- x@upperbound[,k]
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
  if (is.null(main) && length(method <- x@control$name)) {
    #a surveillance algorithm has been run
    action <- switch(class(x), "sts" = "surveillance",
                     "stsNC" = "nowcasting","stsBP" = "backprojection")
    main <- paste0(action, " using ", method)
  }

  #if ylim is not specified, give it a default value
  if(is.null(ylim)){
    ymax <- max(c(observed,upperbound),na.rm=TRUE)
    ylim <- c(-1/20*ymax, ymax)
  }

  # left/right help for constructing the columns
  dx.observed <- 0.5
  upperboundx <- (1:length(upperbound)) - (dx.observed - dx.upperbound)

  #Generate the matrices to plot (values,last value)
  xstuff <- cbind(c(upperboundx,length(observed) + min(1-(dx.observed - dx.upperbound),0.5)))
  ystuff <-cbind(c(upperbound,upperbound[length(observed) ]))

  #Plot the results
  if (length(lty) < 3) lty[3] <- 0  # blank upperbound
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab=ylab,main=main,ylim=ylim,axes = !(xaxis.dates),type=type,lty=lty[-c(1:2)],col=col[-c(1:2)],lwd=lwd[-c(1:2)],...)

  #This draws the polygons containing the number of counts (sep. by NA)
  i <- rep(1:length(observed),each=5)
  dx <- rep(dx.observed * c(-1,-1,1,1,NA), times=length(observed))
  x.points <- i + dx
  y.points <- as.vector(t(cbind(0, observed, observed, 0, NA)))
  polygon(x.points,y.points,col=col[1],border=col[2],lwd=lwd[1],lty=lty[1])

  #Draw upper bound once more in case the polygons are filled
  if (!is.na(col[1])) {
    lines(x=xstuff,y=ystuff,type=type,lty=lty[-c(1:2)],col=col[-c(1:2)],lwd=lwd[-c(1:2)],...)
  }

  #Draw alarm symbols
  alarmIdx <- which(alarm == 1)
  if (length(alarmIdx)>0) {
    matpoints( alarmIdx, rep(-1/40*ylim[2],length(alarmIdx)), pch=alarm.symbol$pch, col=alarm.symbol$col, cex= alarm.symbol$cex, lwd=alarm.symbol$lwd)
  }

  #Draw outbreak symbols
  stateIdx <- which(state == 1)
  if (length(stateIdx)>0) {
    matpoints( stateIdx, rep(-1/20*ylim[2],length(stateIdx)), pch=outbreak.symbol$pch, col=outbreak.symbol$col,cex = outbreak.symbol$cex,lwd=outbreak.symbol$lwd)
  }

  #Label x-axis
  if(xaxis.dates & axes) {
    addFormattedXAxis(x = x, epochsAsDate = epochsAsDate, xaxis.tickFreq = xaxis.tickFreq,
                      xaxis.labelFreq = xaxis.labelFreq, xaxis.labelFormat = xaxis.labelFormat,
                      ...)
  }
  #Label y-axis
  if (axes) {
    axis( side=2 ,...)
  }

  hasupper <- any(upperbound > 0, na.rm = TRUE)
  included <- c(observed = TRUE, upperbound = hasupper,
                outbreaks = length(stateIdx) > 0,
                alarms = hasupper || length(alarmIdx) > 0)
  doLegend <- if (missing(legend.opts)) {
      sum(included) > 1
  } else {
      is.list(legend.opts)
  }
  if(doLegend) {
      ## FIXME: use method-specific default upperbound label?
      ## ublegend <- if (identical(x@control[["ret"]], "value") &&
      ##                 startsWith(x@control$name, "glr"))
      ##                 "GLR statistic" else "Threshold"
      if (!is.na(col[1])) { # filled polygons
          ltyObs <- NA
          pchObs <- 22
          fillObs <- col[1]
      } else {
          ltyObs <- lty[1]
          pchObs <- NA
          fillObs <- NA
      }
      legend.opts <- modifyList(
          list(x = "top",
               lty = c(ltyObs,lty[3],NA,NA)[included],
               lwd = c(lwd[1],lwd[3],outbreak.symbol$lwd,alarm.symbol$lwd)[included],
               col = c(col[2],col[3],outbreak.symbol$col,alarm.symbol$col)[included],
               pch = c(pchObs,NA,    outbreak.symbol$pch,alarm.symbol$pch)[included],
               pt.bg = c(fillObs,NA,NA,NA)[included],
               legend = c("Infected", "Threshold", "Outbreak", "Alarm")[included]),
          legend.opts)
    #Make the legend
    do.call("legend",legend.opts)
  }

  #Call hook function for user customized action using the current environment
  environment(hookFunc) <- environment()
  hookFunc()

  #Extra hook functions for inheritance plotting (see e.g. plot function of stsNC objects)
  environment(.hookFuncInheritance) <- environment()
  .hookFuncInheritance()

  invisible()
}


##############
### alarm plot
##############

stsplot_alarm <- function(
    x, lvl=rep(1,ncol(x)),
    xaxis.tickFreq=list("%Q"=atChange),
    xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="%G\n\n%OQ",
    epochsAsDate=x@epochAsDate, xlab="time", ylab="", main=NULL,
    outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1), # unused
    alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
    cex.yaxis=1, ...)
{
  if (is.null(main) && length(method <- x@control$name)) {
    #a surveillance algorithm has been run
    action <- switch(class(x), "sts" = "surveillance",
                     "stsNC" = "nowcasting","stsBP" = "backprojection")
    main <- paste0(action, " using ", method)
  }

  #Control what axis style is used
  xaxis.dates <- !is.null(xaxis.labelFormat)

  #Initialize plot
  plot(x=c(0.5,nrow(x)+0.5),y=c(0.5,ncol(x)),xlab=xlab,ylab=ylab,main=main,axes=FALSE,type="n",...)

  #Label of x-axis
  if(xaxis.dates){
    addFormattedXAxis(x = x, epochsAsDate = epochsAsDate, xaxis.tickFreq = xaxis.tickFreq,
                      xaxis.labelFreq = xaxis.labelFreq, xaxis.labelFormat = xaxis.labelFormat,
                      ...)
  }
  axis( side=2, at=1:ncol(x),cex.axis=cex.yaxis, labels=colnames(x),las=2)


  #Draw all alarms
  for (i in 1:nrow(x)) {
    idx <- (1:ncol(x))[x@alarm[i,] > 0]
    for (j in idx) {
      points(i,j,pch=alarm.symbol$pch,col=alarm.symbol$col[lvl[j]],cex=alarm.symbol$cex,lwd=alarm.symbol$lwd)
    }
  }

  #Draw lines separating the levels
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
    as.integer(tapply(seq_along(x), INDEX=x, quantile, probs=0.5, type=3))
}
#Only every second unit change
at2ndChange <- function(x,xm1) {
    idxAtChange <- atChange(x,xm1)
    idxAtChange[seq(idxAtChange) %% 2 == 1]
}

#Helper function to format the x-axis of the time series
addFormattedXAxis <- function(x, epochsAsDate = FALSE,
                              xaxis.tickFreq = list("%Q"=atChange),
                              xaxis.labelFreq = xaxis.tickFreq, xaxis.labelFormat = "%G\n\n%OQ",
                              ...) {

  #Old style if there are no Date objects
  if (!epochsAsDate) {
    #Declare commonly used variables.
    nTime <- nrow(x)
    startyear <-  x@start[1]
    firstweek <-  x@start[2]

    if (x@freq ==52) {  #Weekly epochs are the most supported
      # At which indices to put the "at" tick label. This will
      # be exactly those week numbers where the new quarter begins: 1, 14, 27 and 40 + i*52.
      # Note that week number and index is not the same due to the "firstweek" argument
      weeks <- seq_len(nTime) + (firstweek-1)
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
      axis( side=1,labels=FALSE,at=c(1,nTime),lwd.ticks=0,line=1,...)
      axis( at=weekIdx[which(quarter != "I")] , labels=labels.week[which(quarter != "I")] , side=1, line = 1 ,...)
      #Bigger tick marks at the first quarter (i.e. change of the year)
      at <- weekIdx[which(quarter == "I")]
      axis( at=at, labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par("tcl"))

    } else { ##other frequency (not really supported)
      #A label at each unit
      myat.unit <- seq(firstweek,length.out=nTime)

      # get the right year order
      month <- (myat.unit-1) %% x@freq + 1
      year <- (myat.unit - 1) %/% x@freq + startyear
      #construct the computed axis labels -- add quarters if xaxis.units is requested
      mylabels.unit <- paste(year,"\n\n", (myat.unit-1) %% x@freq + 1,sep="")

      #Add axis
      axis( at=seq_len(nTime), labels=NA, side=1, line = 1, ...)
      axis( at=seq_len(nTime)[month==1], labels=mylabels.unit[month==1] , side=1, line = 1 ,...)
      #Bigger tick marks at the first unit
      at <- seq_len(nTime)[(myat.unit - 1) %% x@freq == 0]
      axis( at=at, labels=rep(NA,length(at)), side=1, line = 1 ,tcl=2*par("tcl"))
    }
  } else {
    ################################################################
    #epochAsDate -- experimental functionality to handle ISO 8601
    ################################################################

    dates <- epoch(x, as.Date = TRUE)
    #make one which has one extra element at beginning with same spacing
    datesOneBefore <- c(dates[1]-(dates[2]-dates[1]),dates)

    #Make the line. Use lwd.ticks to get full line but no marks.
    axis( side=1,labels=FALSE,at=c(1,length(dates)),lwd.ticks=0,...)

    ###Make the ticks (depending on the selected level).###
    tcl <- par("tcl")
    tickFactors <- surveillance.options("stsTickFactors")

    #Loop over all pairs in the xaxis.tickFreq list
     for (i in seq_along(xaxis.tickFreq)) {
      format <- names(xaxis.tickFreq)[i]
      xm1x <- as.numeric(formatDate(datesOneBefore,format))
      idx <- xaxis.tickFreq[[i]](x=xm1x[-1],xm1=xm1x[1])

      #Find tick size by table lookup
      tclFactor <- tickFactors[pmatch(format, names(tickFactors))]
      if (is.na(tclFactor)) {
        warning("no \"tcl\" factor found for \"", format ,"\" -> setting it to 1")
        tclFactor <- 1
      }
      axis(1,at=idx, labels=NA,tcl=tclFactor*tcl,...)
    }

    ###Make the labels (depending on the selected level)###
    if (!is.null(xaxis.labelFormat)) {
      labelIdx <- NULL
      for (i in seq_along(xaxis.labelFreq)) {
        format <- names(xaxis.labelFreq)[i]
        xm1x <- as.numeric(formatDate(datesOneBefore,format))
        labelIdx <- c(labelIdx,xaxis.labelFreq[[i]](x=xm1x[-1],xm1=xm1x[1]))
      }

      #Format labels (if any) for the requested subset
      if (length(labelIdx)>0) {
          labels <- rep(NA,nrow(x))
          labels[labelIdx] <- formatDate(dates[labelIdx],xaxis.labelFormat)
          axis(1,at=1:nrow(x), labels=labels,tick=FALSE,...)
      }
    }
  }#end epochAsDate
  #Done
  invisible()
}
