######################################################################
# initialize-method for "stsNC" objects
######################################################################

init.stsNC <- function(.Object, ...,
                       reportingTriangle, predPMF, pi, truth, delayCDF, SR)
{
    .Object <- callNextMethod()  # use initialize,sts-method

    ## initialize defaults for extra stsNC-slots or check supplied values
    dimObserved <- dim(.Object@observed)
    if (missing(pi)) {
        .Object@pi <- array(NA_integer_, dim = c(dimObserved, 2L))
    } else {
        dimPI <- dim(.Object@pi)
        if (length(dimPI) != 3 || any(dimPI != c(dimObserved, 2L)))
            stop("dim(pi) = (", paste0(dimPI, collapse=","), ")")
    }
    if (missing(SR)) {
        .Object@SR <- array(NA_real_, dim = c(nrow(.Object@observed),0L,0L))
    } else {
        stopifnot(length(dim(.Object@SR)) == 3)
    }
    if (missing(truth))
        .Object@truth <- as(.Object, "sts")

    return(.Object)
}

setMethod("initialize", "stsNC", init.stsNC)


######################################################################
# Special coerce method to account for consistent dimensions
######################################################################

setAs(from = "sts", to = "stsNC", function (from) {
    new("stsNC", from,
        pi = array(NA_real_, dim = c(dim(from@observed), 2L)),
        truth = from,
        SR = array(NA_real_, dim = c(nrow(from@observed), 0L, 0L)))
})


######################################################################
# plot-method for the "stsNC" class, which starts by
# using the inherited method, but with some additional plotting
# put into the .hookFunSpecial function.
#
# Parameters:
#  same as the for the plot method of sts objects.
######################################################################

setMethod(f="plot", signature=signature(x="stsNC", y="missing"),
          function (x, type = observed ~ time | unit, ...) {

            ## if special type "delay" (only applies for stsNC objects)
            if (type == "delay") {
              stsNC_plotDelay(x, ...)
              return(invisible())
            }

              ## environment of hook function will be set to evaluation
              ## environment of stsplot_time1() and only then be called
              legend.opts <- lty <- lwd <-
                  "accommodate tools:::.check_code_usage_in_package()"

              #Hook function specifically for nowcasting objects.
              nowcastPlotHook <- function() {
                  #Define some colors for the plotting as well as some plot symbols
                  color <- surveillance.options("colors")
                  pchList   <- c(nowSymbol=10)

                  #Prolong line of last observation (this should go into the plot function
                  idx <- nrow(x) - which.max(!is.na(rev(upperbound(x)))) + 1
                  #Continue line from plot - use same style as stsplot_time1
                  lines( idx+c(-0.5,0.5), rep(upperbound(x)[idx,],2),col=col[3],lwd=lwd[3],lty=lty[3])

                  #Add the prediction intervals as bars (where not NA). Conf level
                  #is found in x@control$alpha
                  idxt <- which(apply(x@pi[1:nrow(x),1,],1,function(x) all(!is.na(x))))
                  for (i in idxt) {
                      lines( i+c(-0.3,0.3), rep(x@pi[i,,1],2),lty=1,col=color["piBars"])
                      lines( i+c(-0.3,0.3), rep(x@pi[i,,2],2),lty=1,col=color["piBars"])
                      lines( rep(i,each=2), x@pi[i,,],lty=2,col=color["piBars"])
                  }

                  #Extract now date and date range of the plotting
                  startDate <- epoch(x)[1]

                  #Add "now" symbol on x-axis. Plotting now takes possible temporal aggregation into account.
                  #points(x@control$now-startDate+1,0,pch=pchList["nowSymbol"],col=color["nowSymbol"],cex=1.5)
                  points(x@control$timeDelay(startDate,x@control$now)+1,0,pch=pchList["nowSymbol"],col=color["nowSymbol"],cex=1.5)
                  #Add this to the legend
                  if (!is.null(legend.opts)) {
                      legend(x="topright",c("Now"),pch=pchList["nowSymbol"],col=color["nowSymbol"],bg="white")
                  }

                  return(invisible())
              }

              callNextMethod(x=x, type=type, ..., .hookFuncInheritance=nowcastPlotHook)
          })



######################################
## For plotting the delay distribution
######################################

######################################################################
## Convert discrete time hazards to PMF
## Parameters:
##  haz - vector with entries for (0,...,Dmax)
######################################################################

haz2pmf <- function(haz) {
    PMF <- 0*haz
    for (i in 0:(length(haz)-1)) {
        PMF[i+1] <- haz[i+1] * (1-sum(PMF[seq(i)]))
    }
    return(PMF)
}

######################################################################
# Find a quantile of a discrete random variable with support on
# 0,...,D and which has a PMF given by the vector prob. We
# define the q quantile as \min_{x} F(x) \geq q.
#
# Parameters:
#   prob - vector on 0,..,D containing the PMF
#   q    - quantile to compute
######################################################################

pmfQuantile <- function(prob,q=0.5) {
    which.max(cumsum(prob) >= q)-1
}


######################################################################
## Show empirical and, if available, model based median of delay
## distribution as a function of occurence time t.
##
## Parameters:
##  nc - nowcast object
##  rT.truth - reporting triangle as it would be at the end. Typically
##             this is taken directly from the nc object.
##  dates - vector of dates where to show the result
##  w - half-width of moving window
##  modelQuantiles - which model quantiles to show
######################################################################

stsNC_plotDelay <- function(nc, rT.truth=NULL, dates=NULL, w=1, modelQuantiles=0.5, epochUnit=NULL) {

  ##Extract reporting triangle from the nc object
  if (is.null(rT.truth)) {
    rT.truth <- reportingTriangle(nc)
  }
  ##Which dates to plot
  if (is.null(dates)) {
    dates <- epoch(nc)
  }

  ##Determine the appropriate unit of the delay
  if (is.null(epochUnit)) {
    epochUnit <- switch( as.character(nc@freq),
                        "12" = "months", "%m" = "months",
                        "52" =  "weeks", "%V"="weeks",
                        "%j"="days", "365" = "days")
  }

  ##Determine max delay from reporting triangle.
  D <- nc@control$D
  res <- matrix(NA, nrow=length(dates), ncol=D+1)

  ##which data variables are actually in rT.truth
  isThere <- !is.na(sapply(dates, function(date) pmatch(as.character(date),rownames(rT.truth))))
  idx <- which(isThere)

  ##Loop over all time points.
  for (i in (w+min(idx)):(max(idx)-w)) {
    now <- dates[i]
    the_idx <- pmatch(as.character(now),rownames(rT.truth))
    subset <- rT.truth[the_idx + c(-w:w),,drop=FALSE]
    res[i,] <- colSums(subset,na.rm=TRUE) / sum(subset,na.rm=TRUE)
  }

  ##A slightly modified function to determine quantiles, which can
  ##handle NAs (if there is no case at all)
  quantile <- function(q) {
    apply(res, 1, function(x) {
      if (all(is.na(x))) return(NA) else return(which.max(cumsum(x) >= q) - 1)
    })
  }

  ##Find 10%, 50% and 90% quantiles
  quants <- sapply(c(0.1,0.5,0.9), quantile)

  ##Make a plot (use plot.Dates instead of matplot)
  plot(dates, quants[,2],xlab="Time of occurence",ylab=paste0("Delay (",epochUnit,")"),ylim=c(0,15),col=1,lty=c(1),lwd=4,type="n")

  idxFirstTruncObs <- which(dates == (nc@control$now - D))
  idxNow <- which(dates == nc@control$now)
  polygon( dates[c(idxFirstTruncObs,idxFirstTruncObs,idxNow,idxNow)], c(-1e99,1e99,1e99,-1e99), col=rgb(0.95,0.95,0.95),lwd=0.001)

  text( dates[round(mean(c(idxNow,idxFirstTruncObs)))], D, "right truncated\n observations",adj=c(0.5,0.5))
  lines(dates, quants[,2],col=1,lty=c(1),lwd=4)
  matlines(dates, quants[,c(1,3)],type="l",col=1,lty=c(2,3),lwd=c(1,1))




  legend_str <- c(expression(q[0.1](T)),expression(q[0.5](T)),expression(q[0.9](T)))
  legend_lty <- c(2,1,3)
  legend_col <- c(1,1,1)
  legend_lwd <- c(1,4,1)
  ##Which dates have been analysed in the nowcasts
  dates2show <- attr(reportingTriangle(nc),"t02s")
  
  ##Loop over all model based estimates
  model_CDF <- delayCDF(nc)
  if (length(model_CDF) > 0) {
    for (methodIdx in seq_len(length(model_CDF))) {
     ##browser()
      ##Fetch CDF from model (can be a vector or a matrix)
      theCDF <- delayCDF(nc)[[names(model_CDF)[methodIdx]]]
      if (!is.matrix(theCDF)) {
        theCDF <- matrix(theCDF, ncol=length(theCDF),nrow=length(dates2show),byrow=TRUE)
      }
      cdf <- cbind(0,theCDF)
      pmf <- t(apply(cdf,1,diff))

      ##Determine model quantiles
      quants.model <- matrix(NA, nrow=length(dates2show),ncol=length(modelQuantiles),dimnames=list(as.character(dates2show),modelQuantiles))
      for (t in 1:length(dates2show)) {
        quants.model[t,] <- sapply(modelQuantiles, function(q) pmfQuantile( pmf[t,],q=q))
      }

      ##Make sure the NAs in the beginning agree
      i <- 1
      while (all(is.na(quants[i,]))) {quants.model[i,] <- NA ; i <- i + 1}


      legend_str <- c(legend_str,substitute(q[0.5]^methodName(T),list(methodName=names(model_CDF)[methodIdx])))
      legend_lty <- c(legend_lty,3+methodIdx)
      legend_col <- c(legend_col,"gray")
      legend_lwd <- c(legend_lwd,2)

      ##only estimates up to 'now' are to be shown and which are within
      ##the moving window of m time points
      show <- (nc@control$now - dates2show <= nc@control$m)
      matlines(dates2show[show], quants.model[show,], col=tail(legend_col,n=1),lwd=ifelse(modelQuantiles==0.5,tail(legend_lwd,n=1),1),lty=ifelse(modelQuantiles==0.5,tail(legend_lty,n=1),2))
    }

    ##Show lines for breakpoints (if available from the model)
    if ("bayes.trunc.ddcp" %in% names(model_CDF)) {
      ddcp.model <- attr(model_CDF[["bayes.trunc.ddcp"]],"model")
      changePoints <- as.Date(colnames(ddcp.model$W))
      for (i in 1:length(changePoints)) {
        axis(1,at=changePoints[i], changePoints[i], las=1, cex.axis=0.7,line=-2.5)
        lines( rep(changePoints[i],2),c(0,1e99),lty=2)
      }
    }
  }

  ##Make a legend
  ##c(expression(q[0.1](T)),expression(q[0.5](T)),expression(q[0.9](T)),expression(q[0.5]^"ddcp"(T)))
  legend(x="bottomleft",legend_str,lty=legend_lty,col=legend_col,lwd=legend_lwd)

  ##Add title
  if (!is.null(nc)) { title(nc@control$now) }

  ##Done
  invisible()
}


