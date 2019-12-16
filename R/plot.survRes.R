plot.survRes.one <- function(x, method=x$control$name, disease=x$control$data, domany=FALSE,ylim=NULL,xaxis.years=TRUE,startyear = 2001, firstweek = 1, xlab="time", ylab="No. infected", main=NULL, type="hhs",lty=c(1,1,2),col=c(1,1,4), outbreak.symbol = list(pch=3, col=3),alarm.symbol=list(pch=24, col=2),legend.opts=list(x="top",legend=c("Infected", "Upperbound", "Alarm", "Outbreak"),lty=NULL,col=NULL,pch=NULL), ...) {

  ################## Handle the NULL arguments ########################################################
  if (is.null(main)) main = paste("Analysis of ", as.character(disease), " using ", as.character(method),sep="")
  #No titles are drawn when more than one is plotted.
  if (domany) main = ""

  survResObj <- x
  observed <- survResObj$disProgObj$observed[survResObj$control$range]
  state    <- survResObj$disProgObj$state[survResObj$control$range]

  #print(list(...))
  # width of the column
  tab <- 0.5

  # left/right help for constructing the columns
  observedxl <- (1:length(observed))-tab
  observedxr <- (1:length(observed))+tab
  upperboundx <- (1:length(survResObj$upperbound)) #-0.5

  # control where the highest value is
  max <- max(max(observed),max(survResObj$upperbound))

  #if ylim is not specified
  #if(is.null(ylim)){
  #  ylim <- c(-1/20*max, max)
  #}

#~~~~~~~~~~~~~~~~~~~~~~~~~~
  if (is.null(ylim)) {
    max <- max(max(observed), max(survResObj$upperbound))
    ylim <- c(-1/20 * max, max)
  } else {
    max <- ylim[2]
  }

  #ensure that there is enough space for the alarm/outbreak symbols
  if(ylim[1]>=0)
    ylim[1] <- -1/20*max
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #Generate the matrices to plot
  xstuff <- cbind(observedxl, observedxr, upperboundx) #no adjusting + min(x$control$range) - 1
  ystuff <- cbind(observed, observed, survResObj$upperbound)


  #Plot the results using one Large plot call (we do this by modifying
  #the call).
  matplot(x=xstuff,y=ystuff,xlab=xlab,ylab=ylab,main=main,ylim=ylim,axes = !(xaxis.years),type=type,lty=lty,col=col,...)


  if (!is.null(survResObj$aggr)) {
    points(upperboundx+tab,survResObj$aggr,col=1)
  }

  for(i in 1:length(observed)){
    matlines( c(i-tab, i+tab), c(observed[i],observed[i]),col=col[1])
    if(survResObj$alarm[i] == 1)
      matpoints( i, -1/40*max, pch=alarm.symbol$pch, col=alarm.symbol$col)
    if(state[i] == 1)
      matpoints( i, -1/20*max, pch=outbreak.symbol$pch, col=outbreak.symbol$col)
  }

  # check where to place the legend. If the left upper side is free place it there
  if (max * 2/3 >= max(
                max(observed[1:floor(1/4 * length(observed))]),
                max(survResObj$upperbound[1:floor(1/4 * length(survResObj$upperbound))])
                )) {
    xlegpos <- 0
  }

  #Label of x-axis
  if(xaxis.years){
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

    # construct the computed axis labels
    #cex <- par()$cex.axis
    #if (cex == 1) {
    mylabels.week <- paste(year,"\n\n",quarter,sep="")
    #} else {
    #  mylabels.week <- paste(year,"\n",quarter,sep="")
    #}

    axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
    axis( side=2 )
  }

  if(is.list(legend.opts)) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$lty)) legend.opts$lty = c(lty[1],lty[3],NA,NA)
    if (is.null(legend.opts$col)) legend.opts$col = c(col[1],col[3],alarm.symbol$col,outbreak.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch = c(NA,NA,alarm.symbol$pch,outbreak.symbol$pch)
    if (is.null(legend.opts$x))   legend.opts$x = "top"
    if (is.null(legend.opts$legend))
      legend.opts$legend = c("Infected", "Upperbound", "Alarm", "Outbreak")

    do.call("legend",legend.opts)
  }


  invisible()
}


#the main function -- cant we do better than this?
plot.survRes <- function(x, method=x$control$name, disease=x$control$data, xaxis.years=TRUE,startyear = 2001, firstweek = 1, same.scale=TRUE,...) {
  observed <- x$disProgObj$observed
  state <- x$disProgObj$state
  alarm <- x$alarm

  #univariate timeseries ?
  if(is.vector(observed))
    observed <- matrix(observed,ncol=1)
  if(is.vector(state))
    state <- matrix(state,ncol=1)
  if(is.vector(alarm))
    alarm <- matrix(alarm,ncol=1)
  nAreas <- ncol(observed)
  max <-  max(max(observed),max(x$upperbound))

  #multivariate time series
  if(nAreas > 1){
    #all areas in one plot
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
        #Create the survRes
        dP <- create.disProg(x$disProgObj$week, observed[,k], state[,k],start=x$start)
        obj <- list(alarm=alarm[,k],disProgObj=dP,control=x$control,upperbound=x$upperbound[,k])
        class(obj) <- "survRes"
        plot.survRes.one(obj,startyear = startyear, firstweek = firstweek,
                         xaxis.years=xaxis.years, ylim=ylim, legend.opts=NULL,domany=TRUE,... )
         mtext(colnames(observed)[k],line=-1.3)
         })
      #reset graphical params
      par(mfrow=c(1,1), mar=c(5, 4, 4, 2)+0.1)
    }
  else {  #univariate time series
    plot.survRes.one(x=x, startyear = startyear, firstweek = firstweek, xaxis.years=xaxis.years, domany=FALSE,...)
  }
  invisible()
}
