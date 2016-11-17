######################################################################
# Function to plot a sequence of nowcasts. Can be wrapped with the
# animation package to produce PDF or Web animations
#
# Parameters:
#  linelist_truth - data.frame containing the linelist of cases/reports
#  dEventCol - name of the column containing the time of event (as Date)
#  dReportCol - name of the column containing the time of report receipt (as Date)
#  aggrgate.by - aggregation level (se function linelist2sts)
#  nowcasts - a list of nowcasts (if NULL then they are generated on the fly - Note: This is currently not implemented!)
#  method - which method to animate. Has to be part of the individual nowcast objects in 'nowcasts'
#  control - control object for controlling how the plotting is done
######################################################################

animate_nowcasts <- function(nowcasts,linelist_truth,
                               method="bayes.trunc.ddcp",
                               control=list(dRange=NULL,anim.dRange=NULL, plot.dRange=NULL,consistent=FALSE,sys.sleep=1,ylim=NULL,cex.names=0.7,col=c("violetred3","#2171B5","orange","blue","black","greenyellow")),showLambda=TRUE) {

  ##Extract the dEventCol and dReportCol from the nowcasts
  dEventCol <- nowcasts[[1]]@control$call$dEventCol
  dReportCol <- nowcasts[[1]]@control$call$dReportCol
  aggregate.by <- nowcasts[[1]]@control$call$aggregate.by

  ##Boolean indicator for those having information on dEventCol
  validVarInfo <- !is.na(linelist_truth[,dEventCol])

  ##Show info about what is being illustrated
  message(paste("Total of ",nrow(linelist_truth)," cases in linelist_truth.\nIllustring reporting for ",sum(!is.na(linelist_truth[,dEventCol]))," cases with information on \"",dEventCol,"\"\n",sep=""))

  ##Reduce linelist_truth to those who have the appropriate information
  linelist_truth <- linelist_truth[validVarInfo,]

  #########################################
  ## Check and set default control options
  #########################################
  if (is.null(control[["dRange",exact=TRUE]])) {
    range <- range(c(linelist_truth[,dEventCol],linelist_truth[,dReportCol]),na.rm=TRUE)
  } else {
    range <- control$dRange
  }
  range.dates <- seq(range[1],range[2],by=aggregate.by)
  #plot.dRange
  if (is.null(control[["plot.dRange",exact=TRUE]])) {
    control$plot.dRange <- range(range)
  }
  #anim.dRange
  if (is.null(control[["anim.dRange",exact=TRUE]])) {
    control$anim.dRange <- control$dRange
  }

  #sys.sleep
  if (is.null(control[["sys.sleep",exact=TRUE]]))
    control$sys.sleep <- 1

  if (is.null(control[["cex.names",exact=TRUE]]))
    control$cex.names <- 1
  if (is.null(control[["col",exact=TRUE]]))
    control$col <- c("violetred3","#2171B5","orange","blue","black","springgreen4")
  if (is.null(control[["showLambda",exact=TRUE]]))
    control$showLambda <- TRUE
  ##Check that a list of nowcasts is available
  if (is.null(nowcasts)) {
    stop("not implemented!")
  }

  #####################
  # Preprocessing block
  #####################

  #Create an sts object with the true number of cases..
  sts <- linelist2sts(linelist_truth,dEventCol,aggregate.by=aggregate.by,dRange=range)

  #Index of the time points in the plot.dRange
  plot.dates.idx <- as.numeric(control$plot.dRange - range[1] + 1)
  #Index of the animate dates
  anim.dates <- seq(control$anim.dRange[1],control$anim.dRange[2],by="1 day")
  idxSet <- pmatch(anim.dates,range.dates)

  ##Find ylim
  if (is.null(control[["ylim",exact=TRUE]])) {
    ymax <- max(observed(sts),upperbound(sts),na.rm=TRUE)
    ymax <- max(ymax,unlist(lapply(nowcasts, function(nc) max(c(observed(nc),upperbound(nc),predint(nc)),na.rm=TRUE))))
    control$ylim <- c(0,ymax)
  }

  ##======================
  ## Loop over all dates
  ##======================

  ##Loop over all days. always show what we know
  for (i in idxSet) { ##fix this
    #Set "today"
    curDate <- as.Date(range.dates[i])
    message("Animating ",as.character(curDate),"...")
    #Choose all reports available until this "today"
    linelist_truth.avail <- linelist_truth[ linelist_truth[,dReportCol] <= curDate,]
    #If consistency checking is requested remove all entries which
    #are "beyond" today
    if (!is.null(control$consistent)) {
      linelist_truth.avail <- linelist_truth.avail[ linelist_truth.avail[,dEventCol] <= curDate,]
    }
    ##Check that date exists in nowcast list.
    sts.nowcast <- nowcasts[[as.character(curDate)]]
    if (is.null(sts.nowcast)) {
      stop("Date: ",as.character(curDate)," not available in nowcasts.")
    }
    ##Check that method exists in nowcast object
    if (!(method %in% nowcasts[[as.character(curDate)]]@control$method)) {
      stop("Method ",method," not in nowcasts!")
    }
    ##Exract the used safePredictLag
    control$safePredictLag <- sts.nowcast@control$now - max(sts.nowcast@control$when)

    ##Fill upperbound and CI slots with output of that method (not pretty code: ToDo Improve!!)
    N.tInf.support <- sts.nowcast@control$N.tInf.support
    Ps <- sts.nowcast@predPMF
    when <- sts.nowcast@control$when
    dateRange <- epoch(sts.nowcast)
    idxt <- which(dateRange %in% when)
    alpha <- sts.nowcast@control$alpha

    ##Loop over all time points
    for (i in idxt) {
      predPMF <- Ps[[as.character(dateRange[i])]][[method]]
      sts.nowcast@upperbound[i,] <- median(N.tInf.support[which.max( cumsum(predPMF)>0.5)])
      sts.nowcast@pi[i,,] <- N.tInf.support[c(which.max(cumsum(predPMF) >= alpha/2),which.max(cumsum(predPMF) >= 1-alpha/2))]
    }
    dimnames(sts.nowcast@pi) <- list(as.character(dateRange),NULL,paste( c(alpha/2*100,(1-alpha/2)*100),"%",sep=""))

    #Done
    upperbound(sts.nowcast)[-idxt] <- NA

    #All events which (in an ideal world) would be available now
    linelist_truth.now <- linelist_truth[ linelist_truth[,dEventCol] <= curDate,]
    sts.now <- linelist2sts(linelist_truth.now,dEventCol,aggregate.by=aggregate.by,dRange=c(range[1],curDate))#range)


    ##Percentage of possible observations which are available
    sum(observed(sts.nowcast))
    sum(upperbound(sts.nowcast))
    message(sprintf("(%.0f%% of total cases in linelist_truth reported)\n",sum(observed(sts.nowcast))/sum(observed(sts.now))*100))

    ##Show the true number of counts
    observed(sts) <- matrix(0,nrow=nrow(sts),ncol=1)
    upperbound(sts) <- matrix(0,nrow=nrow(sts),ncol=1)
    observed(sts)[1:nrow(sts.now),] <- observed(sts.now)
    upperbound(sts)[1:nrow(sts.now),] <- upperbound(sts.now)

    ##Plot the true number of counts as sts object
    plot(sts,legend=NULL,dx.upperbound=0,main="",lwd=c(1,1,3),ylab="No. Cases",ylim=control$ylim,lty=c(1,1,1),axes=FALSE,xlab="",col=c(control$col[c(1,1)],NULL),
         xlim=plot.dates.idx,xaxs="i")

####################start to change. Use proper customizable arguments
### plot the nowcast using the S4 method and then add the other
### stuff on top of it...

    ##Add the nowcast
    plot(sts.nowcast,dx.upperbound=0,axes=FALSE,col=control$col[c(2,2,3)],lty=c(1,1,1),legend=NULL,add=TRUE,lwd=c(3,3,3),xlim=plot.dates.idx,xaxs="i")

    ##Last proper index
    idx <- nrow(sts.nowcast) - which.max(!is.na(rev(upperbound(sts.nowcast)))) + 1
    ##Continue line from plot
    lines( idx+c(-0.5,0.5), rep(upperbound(sts.nowcast)[idx,],2),lty=1,col=control$col[3],lwd=3)

    ##Add CIs from the nowcast
    for (i in 1:nrow(sts.nowcast)) {
      lines( i+c(-0.3,0.3), rep(sts.nowcast@pi[i,,1],2),lty=1,col=control$col[3])
      lines( i+c(-0.3,0.3), rep(sts.nowcast@pi[i,,2],2),lty=1,col=control$col[3])
      lines( rep(i,each=2), sts.nowcast@pi[i,,],lty=2,col=control$col[3])
    }

    ##Add lambda_t if it exists.
    if (method == "bayes.trunc.ddcp" && control$showLambda) {
      lambda <- attr(delayCDF(sts.nowcast)[["bayes.trunc.ddcp"]],"model")$lambda
      showIdx <- seq(ncol(lambda) - control$safePredictLag)
      matlines( showIdx,t(lambda)[showIdx,],col="gray",lwd=c(1,2,1),lty=c(2,1,2))
    }

    ##Add axis information
    axis(2)

    ##Add extra line parts on x-axis
    axis(1,at=0:1e3,tick=TRUE,lwd.ticks=0,labels=rep("",1e3+1))
    axis(1,at=0:1e3,tick=TRUE,lwd.ticks=1,tcl=-0.2,labels=rep("",1e3+1))

    ##Hilight the mondays
    is.monday <- format(range.dates,"%w") == 1
    axis(1,at=(1:length(range.dates))[is.monday],labels=format(range.dates[is.monday],"%a %d %b"),las=2,cex.axis=control$cex.names)

    ##Show month breaks
    dom <- as.numeric(format(range.dates,"%d"))
    axis(1,at=which(dom==1),labels=rep("",sum(dom==1)),tcl=-0.8,lwd=0,lwd.ticks=1)


####################stop to change
    ##Extra
    text <- c("Events up to \"now\"","Reports received by \"now\"",paste("Nowcasts by ",method,sep=""),
              if (method=="bayes.trunc.ddcp") expression(lambda[t]*" of bayes.trunc.ddcp") else NULL)
    col <- c(control$col[1:3],  if (method=="bayes.trunc.ddcp") "gray" else NULL)
    legend(x="topright",text,col=col, lwd=3,lty=1)

    ##Add now symbol
    points(curDate-range[1]+1,0,pch=10,col=control$col[6],cex=1.5)
    ##Add nowcast symbol
    points(curDate-range[1]+1-control$safePredictLag,0,pch=9,col=control$col[3],cex=1.5)
    ##Add this to the legend
    legend(x="right",c("Now","Nowcast horizon"),pch=c(10,9),col=control$col[c(6,3)],pt.cex=1.5)
    ##Pause
    Sys.sleep(control$sys.sleep)
  }

  invisible()
}



