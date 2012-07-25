######################################################################
# Takes a data frame with dates of individual
# cases and create an aggregated sts time series object for these
# data with aggregation occuring at the desired scale.
# Note: This function is experimental
#
# Parameters:
#  linelist - a data frame containing individual case information, one per line
#  dateCol - a character string denoting the column name in case containing
#            the relevant date variable to aggregate
#  aggregate.by - aggregation block length given as a string compatible with
#       seq.Date -- see \link{seq.Date} for further details.
#
# Author: Michael Hoehle
# Date LaMo: 06 Feb 2012
######################################################################

linelist2sts <- function(linelist,dateCol,aggregate.by="1 week",dRange=NULL,
                      startYearFormat=switch(aggregate.by,"1 day"="%V","7 day"="%V","1 week"="%V","1 month"="%Y","3 month"="%Y"),
                      startEpochFormat=switch(aggregate.by,"1 day"="%j","7 day"="%V","1 week"="%V","1 month"="%m","3 month"="%Q")
                      ) {
 #by <- match.arg(by,c("1 day","1 week","1 month","1 year")
  if (is.null(dRange)) {
    dRange <- range(linelist[,dateCol],na.rm=TRUE)
  }

  ## #Make sure that if weeks we span the entire data set.
  ## if ((aggregate.by=="1 week" | aggregate.by == "7 day") & is.null(dRange)) {
  ##   #Adjust first date to a monday and the last to be a sunday
  ##   weekDay <- as.numeric(format(dRange, "%w"))
  ##   dRange[1] <- dRange[1] - ifelse( weekDay[1] == 0, 6, (weekDay[1]-1))
  ##   dRange[2] <- dRange[2] + 7
  ## }

  #Add exactly one time step to dRange to ensure that cut
  #contains the respective level
  maxDate <- seq(max(dRange),length.out=2,by=aggregate.by)[-1]
  dates <- seq(min(dRange), maxDate, by=aggregate.by)

  #Make a table containing the specific number of cases. Note that this
  #needs to occur using a cut statement
  lvl <- cut(linelist[,dateCol], breaks=dates,right=FALSE)

  observed <- table(lvl)
  epoch <- as.Date(names(observed))

  #Translate "by" to freq string
  freq <- switch(aggregate.by,"1 day"=365,"7 day"=52,"1 week"=52,"1 month"=12,"3 month"=4)

  startYear <- as.numeric(formatDate(min(dates),startYearFormat))
  startEpoch <- as.numeric(formatDate(min(dates),startEpochFormat))
                  
  observed <- matrix(observed,ncol=1)

  #Create S4 object
  sts <- new("sts",epoch=as.numeric(epoch),observed=observed, alarm=0*observed, epochAsDate=TRUE,freq=freq,start=c(startYear,startEpoch))

  #Return
  return(sts)
}
