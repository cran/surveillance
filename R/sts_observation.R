################################################################################
### Function for creating an "sts" object with a given observation date
###
### Copyright (C) 2014-2015 Maelle Salmon
################################################################################

sts_observation <- function(sts,dateObservation,cut=TRUE)
{
  # The sts object we shall return
  stsSub <- sts

  # Index of the observation date
  line1 <- which(epoch(sts)==dateObservation)

  # Maximal delay
  D <- dim(stsSub@control$reportingTriangle$n)[2]-1

  # Number of dates
  theEnd <- dim(stsSub@control$reportingTriangle$n)[1]

  # Nothing observed after the observation date (I am a genius)
  stsSub@control$reportingTriangle$n[(line1+1):theEnd,] <- NA
  stsSub@observed[(line1+1):theEnd] <- 0

  # Not everything observed before the observation date

  for (i in 1:D){
    stsSub@control$reportingTriangle$n[line1+1-i,(i+1):(D+1)] <- NA
    stsSub@observed[line1+1-i] <- sum(stsSub@control$reportingTriangle$n[line1+1-i,],na.rm=T)
  }
  stsSub@control$reportingTriangle$n <- stsSub@control$reportingTriangle$n[1:line1,]
  # Return the new sts object
  if (cut){return(stsSub[1:line1])}
  else{return(stsSub)}
}
