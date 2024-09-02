################################################################################
### Simulate count time series with outbreaks (following Noufaily et al, 2012)
###
### Copyright (C) 2014-2015 Maelle Salmon
################################################################################

sts_creation <- function(theta,beta,gamma1,gamma2,m,overdispersion,dates,
                         sizesOutbreak,datesOutbreak,delayMax,alpha,
                         densityDelay)
{
  lengthT <- length(dates)
  # Baseline
  observed <- rep(NA,lengthT)
  upperbound <- rep(NA,lengthT)
  state <- logical(length=lengthT)

  for (t in 1:lengthT) {
    if (m==0){season=0}
    if (m==1){season=gamma1*cos(2*pi*t/52)+ gamma2*sin(2*pi*t/52)}
    if (m==2){season=gamma1*cos(2*pi*t/52)+ gamma2*sin(2*pi*t/52)+gamma1*cos(4*pi*t/52)+ gamma2*sin(4*pi*t/52)}
    mu <- exp(theta + beta*t + season)
    observed[t] <- rnbinom(mu=mu,size=overdispersion,n=1)
    upperbound[t] <- qnbinom(mu=mu,size=overdispersion,p=(1-alpha))

  }

  # Outbreaks
  nOutbreaks <- length(sizesOutbreak)
  if (nOutbreaks>1){
    dens <- lognormDiscrete(Dmax=20,logmu=0,sigma=0.5)
    for (i in 1:nOutbreaks){
      tOutbreak <- which(dates==datesOutbreak[i])
      numberOfCases <- rpois(n=1,lambda=sizesOutbreak[i]*(mu*(1+mu/overdispersion)))
      cases <- rep(0,length(dens))
      if (numberOfCases!=0){
        for (case in 1:numberOfCases){
          t <- sample(x=1:length(dens),size=1,prob=dens)
          cases[t] <- cases[t] + 1
        }
      }
      cases <- cases[cases>0]
      if(sum(cases)>0){
      observed[tOutbreak:(tOutbreak+length(cases)-1)] <- observed[tOutbreak:(tOutbreak+length(cases)-1)] + cases
      state[tOutbreak:(tOutbreak+length(cases)-1)] <- TRUE
      }
    }

  }
  observed <- observed[1:lengthT]

  # Reporting triangle
  if (!is.null(densityDelay)){
    # use density delay
    n <- matrix(0, lengthT, delayMax + 1,dimnames=list(as.character(dates),NULL))
    for (t in 1:lengthT){
      if(observed[t]!=0){
        for (case in 1:observed[t]){
          delay <- sample(x=0:delayMax,size=1,prob=densityDelay)
          if (delay > delayMax) {delay <- delayMax}
          n[t, delay + 1] <- n[t, delay + 1] + 1
        }
      }
    }
  }
  else{
    # Using a poisson as for the outbreaks because it looks good
    n <- matrix(0, lengthT, D + 1,dimnames=list(as.character(dates),NULL))
    for (t in 1:lengthT){
      if(observed[t]!=0){
        for (case in 1:observed[t]){
          delay <- rpois(n=1, lambda=1.5)
          if (delay > D) {delay <- D}
          n[t, delay + 1] <- n[t, delay + 1] + 1
        }
      }
    }
  }
  # Create the sts
  start <- unlist(isoWeekYear(dates[1]), use.names = FALSE)
  newSts <- new("sts", epoch = as.numeric(dates), start = start, upperbound = as.matrix(upperbound),
                freq = 52, observed = observed, state = as.matrix(state), epochAsDate = TRUE)
  newSts@control$reportingTriangle$n <- n
  return(newSts)
}

## FUNCTION FOR DISCRETIZING THE LOG NORM DISTRIBUTION
lognormDiscrete <- function(Dmax=20,logmu=0,sigma=0.5){
  Fd <- plnorm(0:Dmax, meanlog = logmu, sdlog = sigma)
  FdDmax <- plnorm(Dmax, meanlog = logmu, sdlog = sigma)

  #Normalize
  prob <- diff(Fd)/FdDmax
  return(prob)
}
