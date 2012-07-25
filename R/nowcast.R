######################################################################
# Function to perform nowcast at a specific day s. The full
# documentation is available in the nowcast.Rd file.
#
# Parameters:
#  s - a Date object representing today
#  t - a Date object representing the day to do the forecast for.
#      A requirement is that t<=s. ToDo: t should be a vector
#  D - the Database containing columns dHospital and dReport
#
# Returns:
#  stsBP object with estimate and CI in the appropriate slots.
######################################################################

nowcast <- function(s,t,D,dEventCol="dHospital",dReportCol="dReport",
                    method=c("freq.pi","bayes.nb","bayes.betapi","bayes.bnb","uniform"),
                    aggregate.by="1 day",
                    control=list(
                      dRange=NULL,
                      timeDelay=function(d1,d2) {as.numeric(d2-d1)},
                      estimateF="dynamic",
                      alpha=0.05,
                      y.prior.max=300,
                      B=1e5, score=FALSE,PMF=FALSE,sts.truth=FALSE)) {
  
  #Check if t<=s
  if (!all(t<=s)) {
    stop("Assertion t<=s failed.")
  }

  #Check that specified methods are all valid
  method <- match.arg(method,c("freq.pi","bayes.nb","bayes.betapi","bayes.bnb","uniform"),several.ok=TRUE)
  
  
  #If there is a specification of dateRange set dMin and dMax accordingly
  #Otherwise use the limits of the range
  if (is.null(control$dRange)) {
    dMin <- min(D[,dEventCol],na.rm=TRUE)
    dMax <- max(D[,dEventCol],na.rm=TRUE)
  } else {
    dMin <- control$dRange[1]
    dMax <- control$dRange[2]
  }

  dateRange <- seq(dMin,dMax,by=aggregate.by)
  timeDelay <- NULL

  #Choose default timeDelay function if none is given.
  if (is.null(control[["timeDelay"]])) {
    timeDelay <- function(d1,d2) {as.numeric(d2-d1)}
  } else {
    if (is.function(control[["timeDelay"]])) {
      timeDelay <- control[["timeDelay"]]
    } else {
      stop("timeDelay(d1,d2) needs to be a function.")
    }
  }
  
  #Create a column containing the reporting delay using the timeDelay
  #function
  D$delay <- timeDelay(D[,dEventCol],D[,dReportCol])

  #Which observations are available at time s
  D.sub <- D[ na2FALSE(D[,dReportCol] <= s),]
  if (nrow(D.sub)==0) {
    stop(paste("No data available at s=",s,"\n"))
  }
  
  #Create an sts object containing the observed number of counts until s
#  observed <- table(factor(as.character(D.sub[,dEventCol]), levels=as.character(dateRange)))
#  sts <- new("sts",epoch=as.numeric(dateRange),observed=matrix(observed,ncol=1),epochAsDate=TRUE,freq=365)
  sts <- linelist2sts(D.sub,dEventCol,aggregate.by=aggregate.by,dRange=dateRange)
  sts <- as(sts,"stsBP")

  #Create an object containing the "truth" based on D
#  observed <- table(factor(as.character(D[,dEventCol]), levels=as.character(dateRange)))
#  sts.truth <- new("sts",epoch=as.numeric(dateRange),observed=matrix(observed,ncol=1),epochAsDate=TRUE,freq=365)
  sts.truth <- linelist2sts(D,dEventCol,aggregate.by=aggregate.by,dRange=dateRange)
  if (is.null(control[["sts.truth",exact=TRUE]])) {
    control$sts.truth <- FALSE
  }

  #Estimation function for the delay. Standard procedure is to reduce
  #database to only contain the cases which are available at time s
  if (is.null(control[["estimateF",exact=TRUE]]) || (control[["estimateF",exact=TRUE]] == "dynamic")) {
    estimateF <- function(Ds,s,dReportCol) {
      if (nrow(Ds)>0) {
        F <- ecdf(Ds$delay)
      } else {
        warning("No data available to estimate the ECDF at time ",as.character(s),"\n")
        F <- stepfun(x=0,y=c(0,0))
      }
      attr(F,"ms") <- sum(!is.na(Ds$delay))
      return(F)
    }
  } else {
    if (!is.function(control[["estimateF",exact=TRUE]])) {
      stop("The argument control$estimateF needs to be a function.")
    }
  }
  #Estimate delay CDF at s using the specified function
  F <- estimateF(D.sub,s,dReportCol)
  #How many observations were used for the estimation
  ms <- attr(F,"ms")
  
  ######################################################################
  #y.prior max including check that y.prior.max is not too small.
  ######################################################################
  if (is.null(control[["y.prior.max",exact=TRUE]])) {
    control$y.prior.max=300
  }
  y.prior.max <- control$y.prior.max
  if (2*y.prior.max < max(observed(sts),na.rm=TRUE)) {
    warning("y.prior.max appears too small. Largest observed value is more than 50% of y.prior.max, which -- in case this number is extrapolated -- might cause problems.\n")
  }
  #Create a vector representing the support of y_{t}
  yt.support <- 0:y.prior.max
  
  if (is.null(control[["alpha",exact=TRUE]])) {
    alpha <- 0.05
  } else {
    alpha <- control$alpha
  }
  if (is.null(control[["B",exact=TRUE]])) {
    B <- 1e5
  } else {
    B <- control$B
  }
  if (is.null(control[["score",exact=TRUE]])) {
    control$score <- FALSE
  } 

  #List of scores to calculate. Can become an argument later on
  scores <- c("logS","RPS","dist.median","outside.ci")

  #Initialize scoring rule results - to be saved in control slot -- dirty
  SR <- array(0,dim=c(nrow(sts),length(method),length(scores)))

  #List for storing the PMFs.
  if (is.null(control[["PMF",exact=TRUE]])) {
    control$PMF <- FALSE
  }

  ######################################################################
  # Done manipulating the control list with default arguments
  ######################################################################
  sts@control <- control
  
  #Save truth 
  if (control$sts.truth) {
    sts@control$sts.truth <- sts.truth
  } else {
    sts@control$sts.truth <- NULL
  }

  sts@control$PMF <- list()
  
  ######################################################################
  # Helper functions for Bayesian now-casting.
  ######################################################################

  #PMF of the Discrete uniform between lower and upper
  ddu <- function(x,lower=0,upper=1) {
    if (lower>upper) stop("lower > upper")
    return( 1/(upper-lower+1) * ((x>= lower) & (x<= upper)))
  }
    
  #Beta distribution after having observed a sample of x_s positives and
  #m_s - x_s negatives. Protect especially against 0/1 evaluation
  dpits <- function(pits) {
    val <- dbeta(pits, shape1=1/2 + F(diffsti)*ms, shape2=1/2 + ms -  F(diffsti)*ms)
    ifelse(is.finite(val),val,1e99)
  }

  ######################################################################
  # Posterior based on the beta-negative binomial distribution
  ######################################################################

  dbnb <- function(k,n,alpha,beta) {
    #Check if k's outside the support are requested.
    neg <- k<0
    k[neg] <- 0
    #Calculate the density of the beta-negbin. See Teerapabolarn (2008)
    num <- lgamma(n+alpha)+lgamma(k+beta)+lgamma(n+k)+lgamma(alpha+beta)
    den <- lgamma(n+k+alpha+beta)+lgamma(n)+lgamma(k+1)+lgamma(alpha)+lgamma(beta)
    res <- exp(num-den)
    res[neg] <- 0
    return( res)
  }
  
  dpost.bnb <- function(yt) {
    #Shape parameters of the beta distribution for the proportion p
    alpha <- 1/2 + F(diffsti)*ms
    beta  <- 1/2 + ms -  F(diffsti)*ms
    #Add yt, because the support above is 0,... (i.e. only failures are counted)
    return(dbnb( yt-yts,n=yts+1,alpha=alpha,beta=beta))
  }
  
  ######################################################################
  # Posterior based on negative binomial approximation
  #  yt - vector of integer values where to evaluate the function
  #  pits - the probability
  ######################################################################
  dpost.nb <- function(yt,pits) {
    ifelse(yt < yts,0, dnbinom( yt-yts, yts+1 , prob=pits))
  }
  
  #Unormalized joint posterior of N and \pi_{ts}
  post.ytpi.unorm <- function(yt,yts,pits) {
    dbinom( yts, size=yt, prob=pits) * ddu(yt,upper=y.prior.max) * dpits(pits)
  }
  
      
  ######################################################################
  #Loop over all time points in t
  ######################################################################
  idxt <- which(dateRange %in% t)
  for (i in idxt) {

    #Calculate time difference between dateRange[i] and s, i.e. typically
    #s-dateRange[i] once and for all
    diffsti <- timeDelay(dateRange[i], s)

    #Observed value
    yts <- as.numeric(observed(sts)[i,])

    #Calculate pi estimates (fixed & dynamic)
    pits <- F(diffsti)

    #List of casts containing probability distributions
    Ps <- list()
    
    ######################################################################
    # Simple frequentist now-casting using pointwise logit bounds
    ######################################################################

    if ("freq.pi" %in% method) {
      #Proportion reported
      pi.hat <- F(diffsti)# mean(with(D.sub, delay <= diffsti),na.rm=TRUE)
      
      #Use that on logit scale we have an asymptotic normal distribution.
      #See, e.g., Lachin (2000), p.17. -- here dynamic estimation
      if (pi.hat != 0 & pi.hat != 1) { 
        se.logit.pi.hat <- sqrt(1/(ms*pi.hat*(1-pi.hat)))
        #"Frequentisish" predictive distribution approach: based on
        #asymptotic normal on logit transform.
        #This only takes uncertainty of estimation into account
        pits.sample <- plogis(rnorm(B , qlogis(pi.hat), sd=se.logit.pi.hat))
        yts.sample <- round(observed(sts)[i,]/pits.sample)
      } else { #no uncertainty
        yts.sample <- rep(observed(sts)[i,],B)
      }
      #PMF of the frequentist approach. Add 1/B to all cells to make sure there are no zeroes (HACK!)
#      (table(factor(yts.sample, levels=yt.support))+1)/(B+length(yt.support))
      Ps[["freq.pi"]] <- table(factor(yts.sample, levels=yt.support))/B
    }

    ######################################################################
    #For Negbin approximation the posterior is available, but for consistency
    #we create a vector with the values of the PMF
    ######################################################################

    if ("bayes.nb" %in% method) {
      Ps[["bayes.nb"]] <- dpost.nb(yt.support,pits)
    }
  
    ######################################################################
    #A really bad forecast -- the uniform
    ######################################################################
    if ("uniform" %in% method) {
      Ps[["uniform"]] <- rep(1/length(yt.support),length(yt.support))
    }
    
    ######################################################################
    #Posterior which does not take the uncertainty in the estimation
    #of \pi_{ts} into account. This can be approximated by a NegBin
    ######################################################################

#    P.man.dynamic.pifix <- post.ytpi.unorm(yt=yt.support, pits=pits)
#    P.man.dynamic.pifix <- P.man.dynamic.pifix/sum(P.man.dynamic.pifix)

    ######################################################################
    #Alternative. Incorporate uncertainty of the estimation of pits
    ######################################################################
    if ("bayes.betapi" %in% method) {
      P.bayes.betapi <- numeric(length(yt.support))
      for ( j in 1:length(yt.support)) {
        yt <- yt.support[j]
        f <- function(pits) { post.ytpi.unorm(yt=yt, yts=yts, pits=pits) }
        P.bayes.betapi[j] <- integrate(f,lower=0,upper=1)$value
      }
      Ps[["bayes.betapi"]] <- P.bayes.betapi/sum(P.bayes.betapi)
    }

    ######################################################################
    # Beta-negative binomial. Should be the same as the bayes.betapi
    ######################################################################

    if ("bayes.bnb" %in% method) {
      Ps[["bayes.bnb"]] <- dpost.bnb(yt.support)
    }

    #Compute as suggested in Brookmeyer and Damiano (1989)
    #if ("brookmeyer.damiano" %in% method) {
    #  
    #}
    
    ######################################################################
    ######################################################################
    #Done with the computation of PMFs. Now use these for calculation.
    ######################################################################
    ######################################################################
    
    #Save PMFs if thats requested
    if (control$PMF) {
      sts@control$PMF[[as.character(t)]] <- Ps
    }
    
    #Evaluate scoring rules, if requested 
    if (control$score) {
      #Infer the true value
      ytinf <- observed(sts.truth)[i,]

      #Evaluate all scores for all predictive distributions
      for (i.P in 1:length(Ps)) {
        for (i.score in 1:length(scores)) {
          SR[i,i.P,i.score] <- do.call(scores[i.score],args=list(P=Ps[[i.P]],y=ytinf,alpha=alpha))
        }
      }
    } #end if control$score

    #Add first cast & ci to stsBP slots
    sts@upperbound[i,] <- median(yt.support[which.max( cumsum(Ps[[1]])>0.5)])
    sts@ci[i,,] <- yt.support[c(which.max(cumsum(Ps[[1]]) > alpha/2),which.max(cumsum(Ps[[1]]) > 1-alpha/2))]
    

  } #end of loop over time points

  #Add scoring rule to output
  if (control$score) {
    dimnames(SR)    <- list(as.character(dateRange),names(Ps),scores)
    sts@control$SR <- SR
  } else {
    sts@control$SR <- NULL
  }
  #Other arguments to save
  sts@control$yt.support <- yt.support
  
  #Done
  return(sts)
}

## #Example section
## example <- function() {
##   library("surveillance")

##   #Get some data from somewhere
##   source("nowcast2.R")
##   D <- loadData()

##   #Test the function
##   source("nowcast-surveillance.R")
##   s <- as.Date("2011-06-02") ;
##   k <- 10
##   l <- 3
##   t <- seq(s-k-l+1,s-l,by="1 day")
##   dRange <- as.Date(c("2011-05-01","2011-07-10"))
##   nc1 <- nowcast(s=s,t=t,D=D,method="bayes.nb",control=list(dRange=dRange,score=TRUE))

##   #Sow result
##   plot(nc1,xaxis.years=FALSE,dx.upperbound=0,legend=NULL,lty=c(1,1,1),lwd=c(1,1,2),ylab="Cases",xlab="Time (days)",main="")
##   idx <- max(which(!is.na(upperbound(nc1))))
##   lines( c(idx-0.5,idx+0.5), rep(upperbound(nc1)[idx,],2),lwd=2,col="blue")
  
##   ##Show CIs
##   for (i in 1:nrow(nc1)) {
##     points(i, upperbound(nc1)[i,], col="indianred")
##     lines( i+c(-0.3,0.3), rep(nc1@ci[i,,1],2),lty=1,col="indianred2")
##     lines( i+c(-0.3,0.3), rep(nc1@ci[i,,2],2),lty=1,col="indianred2")
##     lines( rep(i,each=2), nc1@ci[i,,],lty=2,col="indianred2")
##   }
##   #Add "now" on the x-axis
##   points( as.numeric(s-dRange[1])+1,0,pch=10,cex=1.5,col="red")


##   #Same as animation
## #  scoreRange <- seq(as.Date("2011-05-25"),max(dRange),by="1 day")
##   scoreRange <- seq(as.Date("2011-05-15"),max(dRange),by="1 day")
##   for (i in 1:length(scoreRange)) {
##     s <- scoreRange[i]
##     t <- seq(s-k-l+1, s-l, by="1 day")
##     nc1 <- nowcast(s=s,t=t,D=D,method="bayes.nb",control=list(dRange=dRange))

##     #Sow result
##     plot(nc1,xaxis.years=FALSE,dx.upperbound=0,legend=NULL,lty=c(1,1,1),lwd=c(1,1,2),ylab="Cases",xlab="Time (days)",main="",ylim=c(0,80))
##     idx <- max(which(!is.na(upperbound(nc1))))
##     lines( c(idx-0.5,idx+0.5), rep(upperbound(nc1)[idx,],2),lwd=2,col="blue")
  
##      ##Show CIs
##     for (i in 1:nrow(nc1)) {
##       points(i, upperbound(nc1)[i,], col="indianred")
##       lines( i+c(-0.3,0.3), rep(nc1@ci[i,,1],2),lty=1,col="indianred2")
##       lines( i+c(-0.3,0.3), rep(nc1@ci[i,,2],2),lty=1,col="indianred2")
##       lines( rep(i,each=2), nc1@ci[i,,],lty=2,col="indianred2")
##     }
##     #Add "now" on the x-axis
##     points( as.numeric(s-dRange[1])+1,0,pch=10,cex=1.5,col="red")
##     Sys.sleep(0.5)
##   }
## }


######################################################################
# Helper functions
######################################################################

#Helper function
na2FALSE <- function(x) {x[is.na(x)] <- FALSE ; return(x) }

######################################################################
# Logarithmic score
#
# Parameters:
#  P - predictive distribution, given as a vector containing the PMF
#      with support 0,...,N.prior.max
#  y - the actual observation. Can be a vector.
#
# Returns:
#  -log P(y). If y outside 0,..,N.prior.max then -Inf.
######################################################################

logS <- function(P, y, ...) {
  return(ifelse( y>=0 & y<=length(P)-1, -log(P[y+1]), -Inf))
}

######################################################################
# Ranked probability score
#
# Parameters:
#  P - predictive distribution, given as a vector containing the PMF
#      with support 0,...,N.prior.max
#  y - the actual observation. Can be a vector.
#
# Returns:
#  -log P(y). If y outside 0,..,N.prior.max then -Inf.
######################################################################

RPS <- function(P,y, ...) {
  N.support <- 0:(length(P)-1)
  sum( (cumsum(P) -  (y <= N.support))^2)
}

#Some other scoring rules which are not proper.
dist.median <- function(P,y, ...) {
  point.estimate <- which.max(cumsum(P)>0.5) - 1
  return(abs(point.estimate - y))
}

#0/1 indicator of observed value outside equal tailed (1-alpha/2) CI
outside.ci <- function(P,y,alpha) {
  N.support <- 0:(length(P)-1)
  ci <- N.support[c(which.max(cumsum(P) > alpha/2),which.max(cumsum(P) >
1-alpha/2))]
  ifelse( y>=ci[1] & y<=ci[2], 0, 1)
}
