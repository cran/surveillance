###################################################
### chunk number 1:
###################################################
create.disProg <- function(week, observed, state, start=c(2001,1), freq=52, neighbourhood=NULL, populationFrac=NULL,epochAsDate=FALSE){
  namesObs <-colnames(observed)

  # check whether observed contains only numbers
  if(!all(sapply(observed, is.numeric))){
    stop("\'observed\' must be a matrix with numbers\n")
  }

  #univariate timeseries ?
  if(is.vector(observed)){
    observed <- matrix(observed,ncol=1)
    namesObs <- deparse(quote(observed))
  } else {  # ensure we have a matrix
    observed <- as.matrix(observed)
  }

  if(missing(state)){
    state <- 0*observed
  } else if(is.vector(state)){
    state <- matrix(state,ncol=1)
  } else {
    state <- as.matrix(state)
  }

  #check number of columns of observed and state
  nAreas <- ncol(observed)
  nObs <- nrow(observed)
  if(ncol(observed) != ncol(state)){
    #if there is only one state-vector for more than one area, repeat it
    if(ncol(state)==1) {
      state <- matrix(rep(state,nAreas),ncol=nAreas,byrow=FALSE)
    } else {
      cat('wrong dimensions of observed and state \n')
      return(NULL)
    }
  }

  #check neighbourhood matrix
  # neighbourhood can be a matrix or an array of dimension c(nAreas,nAreas, nrow(observed))
  if(!is.null(neighbourhood) ) {
    dimNhood <- dim(neighbourhood)
    if(length(dimNhood)==2 & any(dimNhood != nAreas)) {
      cat('wrong dimensions of neighbourhood matrix \n')
      return(NULL)
    } else if (length(dimNhood)==3 & (any(dimNhood[1:2] != nAreas) | (dimNhood[3] != nrow(observed)) )){
      cat('wrong dimensions of neighbourhood matrix \n')
      return(NULL)
    }
  } else {
     # no neighbourhood specified
     neighbourhood <- matrix(NA,nrow=nAreas,ncol=nAreas)
  }

  if(is.null(populationFrac)) {
    populationFrac <- matrix(1/ncol(observed),nrow=nObs, ncol=ncol(observed))
  } else {
    # make sure populationFrac is a matrix
    populationFrac <- as.matrix(populationFrac)
    # check dimensions
    if(nrow(populationFrac)!= nObs | ncol(populationFrac)!= nAreas)
      stop("dimensions of \'populationFrac\' and \'observed\' do not match\n")
    # check whether populationFrac contains only numbers
    if(!all(sapply(populationFrac, is.numeric))){
      stop("\'populationFrac\' must be a matrix with real numbers\n")
  }
  }

  #labels for observed and state
  if(is.null(namesObs)){
    namesObs <- paste(deparse(quote(observed)),1:nAreas,sep="")

  }

  colnames(observed) <- namesObs
  colnames(state) <- namesObs

  res <- list("week"=week, "observed"=observed, "state"=state, "start"=start, "freq"=freq,  "neighbourhood"=neighbourhood, "populationFrac"=populationFrac,"epochAsDate"=epochAsDate)
  class(res) <- "disProg"
  return(res)
}

print.disProg <- function(x, ...) {
  cat( "-- An object of class disProg -- \n" )
  cat( "freq:\t\t", x$freq,"\n" )
  cat( "start:\t\t", x$start,"\n" )
  cat( "dim(observed):\t", dim(x$observed), "\n\n")

  n <- 1
  cat("Head of observed:\n")
  print(head(x$observed,n))

  #cat("\nhead of neighbourhood:\n")
  #print( head(x$neighbourhood,n))
}


###################################################
### chunk number 3:
###################################################
aggregate.disProg <- function(x,...){
  #aggregate observed counts
  observed <- apply(x$observed,MARGIN=1,sum)
  #aggregate states
  state <- apply(x$state,MARGIN=1,sum)
  state[state > 1] <- 1

  #create univariate disProg object
  x <- create.disProg(week=x$week, observed=observed, state=state, freq=x$freq,start=x$start)
  return(x)
}


###################################################
### chunk number 4:
###################################################

plot.disProg.one <- function(x, title = "", xaxis.years=TRUE, quarters=TRUE, startyear = x$start[1], firstweek = x$start[2], ylim=NULL, xlab="time", ylab="No. infected",type="hh",lty=c(1,1),col=c(1,1), outbreak.symbol = list(pch=3, col=3),legend.opts=list(x="top", legend=c("Infected", "Outbreak"),lty=NULL,pch=NULL,col=NULL),...) {


  observed <- x$observed
  state    <- x$state

  # width of the column
  tab <- 0.5

  # left/right help for constructing the columns
  observedxl <- (1:length(observed))-tab
  observedxr <- (1:length(observed))+tab

  # control where the highest value is
  max <- max(observed)

  #if ylim is not specified
  if(is.null(ylim)){
    ylim <- c(-1/20*max, max)
  }

  #Plot the results using one Large plot call
  matplot(x=cbind(observedxl, observedxr),y=cbind(observed, observed),xlab=xlab,ylab=ylab,
          type=type,lty=lty, col=col, ylim=ylim,axes = !(xaxis.years),...)

  #Show the outbreaks
  if (!is.null(outbreak.symbol)) {
    for(i in 1:length(observed)){
      matlines( c(i-tab, i+tab), c(observed[i],observed[i]) )
      if(state[i] == 1)
        matpoints( i, ylim[1], pch=outbreak.symbol$pch, col=outbreak.symbol$col)
    }
  }

  title(title)
  cex <- par()$cex.axis

  #Label of x-axis
  if(xaxis.years){
    # get the number of quarters lying in range for getting the year and quarter order
    obsPerYear <- x$freq
    obsPerQuarter <- x$freq/4
    myat.week <- seq(ceiling((obsPerYear-firstweek+1)/obsPerQuarter) * obsPerQuarter + 1, length(observed)+(floor((obsPerYear-firstweek + 1)/obsPerQuarter) * obsPerQuarter +1), by=obsPerQuarter)
    # get the right year order
    year <- (myat.week - obsPerYear) %/% obsPerYear + startyear
    # function to define the quarter order
    quarterFunc <- function(i) { switch(i+1,"I","II","III","IV")}
    # get the right number and order of quarter labels
    quarter <- sapply( (myat.week-1) %/% obsPerQuarter %% 4, quarterFunc)
    # get the positions for the axis labels
    myat.week <- myat.week - (obsPerYear - firstweek + 1)

    # construct the computed axis labels
    if (quarters) {
      if (cex == 1) {
        mylabels.week <- paste(year,"\n\n",quarter,sep="")
      } else {
        mylabels.week <- paste(year,"\n",quarter,sep="")
      }
    } else {
      mylabels.week <- paste(year,sep="")
    }


    axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
    axis( side=2 )
  }

  #should there be a legend?
  if(is.list(legend.opts)) {
    #Fill empty (mandatory) slots in legend.opts list
    if (is.null(legend.opts$lty)) legend.opts$lty = c(lty[1],NA)
    if (is.null(legend.opts$col)) legend.opts$col = c(col[1],outbreak.symbol$col)
    if (is.null(legend.opts$pch)) legend.opts$pch = c(NA,outbreak.symbol$pch)
    if (is.null(legend.opts$x))   legend.opts$x = "top"
    if (is.null(legend.opts$legend)) legend.opts$legend = c("Infected", "Outbreak")

    #Create the legend
    do.call("legend",legend.opts)
  }

  invisible()
}

plot.disProg <- function(x, title = "", xaxis.years=TRUE, startyear = x$start[1], firstweek = x$start[2], as.one=TRUE, same.scale=TRUE, ...){
  if (xaxis.years && isTRUE(x[["epochAsDate"]]))
    warning("plot.disProg can't handle Date entries; axis labels are based on 'start'")

  observed <- x$observed
  state    <- x$state

  #univariate timeseries ?
  if(is.vector(observed))
    observed <- matrix(observed,ncol=1)
  if(is.vector(state))
    state <- matrix(state,ncol=1)

  nAreas <- ncol(observed)
  max <- max(observed)

  #check if x is multivariate or univariate

  #multivariate time series
  if(nAreas > 1){
    #all areas in one plot -- not supported in sts
    if(as.one){
      matplot(observed,type="l",lty=1:nAreas,col=1:nAreas,ylim=c(0, 1.1*max),xlab="time",ylab="No. of Infected", axes=!xaxis.years)
      #If no legend.opts is specified or not set to null
      if ((is.na(pmatch("legend.opts",names(list(...))))) |
          (!is.na(pmatch("legend.opts",names(list(...)))) & (!is.null(list(...)$legend.opts)))) {
        legend.opts <- list(...)$legend.opts
        if (is.null(legend.opts$x)) legend.opts$x = "topleft"
        if (is.null(legend.opts$legend)) legend.opts$legend = colnames(observed)
        if (is.null(legend.opts$col)) legend.opts$col = 1:nAreas
        if (is.null(legend.opts$lty)) legend.opts$lty = 1:nAreas
        if (is.null(legend.opts$ncol)) legend.opts$ncol = 5
        if (is.null(legend.opts$bty)) legend.opts$bty = "n"

        do.call("legend",legend.opts)
      }

      title(title)

      if(xaxis.years){  #todo: move this as output of ONE function
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
          cex <- par()$cex.axis
          if (cex == 1) {
            mylabels.week <- paste(year,"\n\n",quarter,sep="")
          } else {
            mylabels.week <- paste(year,"\n",quarter,sep="")
          }

          axis( at=myat.week , labels=mylabels.week , side=1, line = 1 )
          axis( side=2 )
       }

    } else {  #plot each area
      #set window size
      par(mfrow=magic.dim(nAreas),mar=c(2,1,1,1))

      if(same.scale)
        ylim <- c(-1/20*max, max)
      else
        ylim <- NULL

      #plot areas
      k <- 1:nAreas
      sapply(k, function(k) {
         plot.disProg.one(create.disProg(x$week, observed[,k], state[,k], freq=x$freq,start=x$start),
                          title = "", startyear = startyear, firstweek = firstweek,
                          xaxis.years=xaxis.years, ylim=ylim, legend.opts=NULL, ... )
         mtext(colnames(observed)[k],line=-1.3)
         })
      #reset graphical params
      par(mfrow=c(1,1), mar=c(5, 4, 4, 2)+0.1)
    }
  } else {  #univariate time series
    plot.disProg.one(x=x, title = title, startyear = startyear, firstweek = firstweek, xaxis.years=xaxis.years, ...)
  }
  invisible()
}
