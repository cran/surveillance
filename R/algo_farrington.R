###################################################
### chunk number 1: 
###################################################
anscombe.residuals <- function(m,phi) {
  y <- m$y
  mu <- fitted.values(m)
  
  #Compute raw Anscombe residuals
  a <- 3/2*(y^(2/3) * mu^(-1/6) - mu^(1/2))
  #Compute standardized residuals
  a <- a/sqrt(phi * (1-hatvalues(m)))
  return(a)
}


###################################################
### chunk number 2: 
###################################################
algo.farrington.assign.weights <- function(s) {
  #s_i^(-2) for s_i<1 and 1 otherwise
  gamma <- length(s)/(sum(  (s^(-2))^(s>1) ))
  omega <- numeric(length(s)) 
  omega[s>1] <- gamma*(s[s>1]^(-2))
  omega[s<=1] <- gamma
  return(omega)
}


###################################################
### chunk number 3: 
###################################################
algo.farrington.fitGLM <- function(response,wtime,timeTrend=TRUE,reweight=TRUE) {
  #Model formula depends on whether to include a time trend or not.
  theModel <- as.formula(ifelse(timeTrend, "response~1+wtime","response~1"))

  #Fit it.
  model <- glm(theModel, family = quasipoisson(link="log"))
    
 #Check convergence - if no convergence we return empty handed.
  if (!model$converged) {
    #Try without time dependence
    if (timeTrend) {
     model <- glm(response ~ 1, family = quasipoisson(link="log"))
     cat("Warning: No convergence with timeTrend -- trying without.\n")
    } 

    if (!model$converged) {
      cat("Warning: No convergence in this case.\n")
      print(cbind(response,wtime))
      return(NULL)
    }
  }

  #Overdispersion parameter phi
  phi <- max(summary(model)$dispersion,1)
  
  #In case reweighting using Anscome residuals is requested
  if (reweight) {
    s <- anscombe.residuals(model,phi)
    omega <- algo.farrington.assign.weights(s)
    model <- glm(theModel,family=quasipoisson(link="log"),weights=omega)
    #Here, the overdispersion often becomes small, so we use the max
    #to ensure we don't operate with quantities less than 1.
    phi <- max(summary(model)$dispersion,1)
  } # end of refit.
  

  #Add wtime, response and phi to the model
  model$phi <- phi
  model$wtime <- wtime
  model$response <- response
  #Done
  return(model)
}


###################################################
### chunk number 4: 
###################################################

algo.farrington.threshold <- function(pred,phi,alpha=0.01,skewness.transform="none") {
  #Fetch mu0 and var(mu0) from the prediction object
  mu0 <- pred$fit
  tau <- phi + (pred$se.fit^2)/mu0
  #Standard deviation of prediction, i.e. sqrt(var(h(Y_0)-h(\mu_0))) 
  switch(skewness.transform,
         "none" = { se <- sqrt(mu0*tau); exponent <- 1},
         "sqrt" = { se <- sqrt(1/4*tau); exponent <- 1/2},
         "2/3"  = { se <- sqrt(4/9*mu0^(1/3)*tau); exponent <- 2/3})

  #Note that lu can contain NA's if e.g. (-1.47)^(3/2)
  lu <- sort((mu0^exponent + c(-1,1)*qnorm(1-alpha/2)*se)^(1/exponent),na.last=FALSE)

  #Ensure that lower bound is non-negative
  lu[1] <- max(0,lu[1],na.rm=TRUE)

  #Return lower and upper bounds
  return(lu)
}


###################################################
### chunk number 5: 
###################################################
algo.farrington <- function(disProgObj, control=list(range=NULL, b=3, w=3, reweight=TRUE, verbose=FALSE,alpha=0.01)) { 
  #Fetch observed
  observed <- disProgObj$observed
  freq <- disProgObj$freq

  ######################################################################
  # Fix missing control options
  ######################################################################
  if (is.null(control$range)) {
    control$range <- (freq*control$b - control$w):length(observed)
  }
  if (is.null(control$b))        {control$b=5}
  if (is.null(control$w))        {control$w=3}
  if (is.null(control$reweight)) {control$reweight=TRUE}
  if (is.null(control$verbose))  {control$verbose=FALSE}
  if (is.null(control$alpha))    {control$alpha=0.05}
  if (is.null(control$trend))    {control$trend=TRUE}
  if (is.null(control$plot))     {control$plot=FALSE}
  if (is.null(control$limit54))  {control$limit54=c(5,4)}

  #check options
  if (!((control$limit54[1] >= 0) &  (control$limit54[2] > 0))) {
    stop("The limit54 arguments are out of bounds: cases >= 0 and perior > 0.")
  }

  # initialize the necessary vectors
  alarm <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  trend <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  upperbound <- matrix(data = 0, nrow = length(control$range), ncol = 1)
  
  # Define objects
  n <- control$b*(2*control$w+1)
  
  # 2: Fit of the initial model and first estimation of mean and dispersion
  #    parameter
  for(k in control$range) {
    # transform the observed vector in the way
    # that the timepoint to be evaluated is at last position
    #shortObserved <- observed[1:(maxRange - k + 1)]

    if (control$verbose) { cat("k=",k,"\n")}

    #Find all weeks with index k-w,..,k+w 
    #in the years (current year)-1,...,(current year)-b
    wtime <- NULL
    for (i in control$b:1){
      wtime <- append(wtime,seq(k-freq*i-control$w,k-freq*i+control$w,by=1))
    }
    response <- NULL # die Responsespalte
    for (i in (control$b:1)) {
      if (control$verbose) {cat("b=",i,"\trange=",((k-i*freq)-control$w):((k-i*freq)+control$w),"\n")}

      for (j in (((k-i*freq)-control$w):((k-i*freq)+control$w))){
        if (j<1) {
          cat("Warning: Selection index less than 1!\n")
        }
        else {
          response <- append(response,observed[j])
        }
      }
    } 
    if (control$verbose) { print(response)}

    ######################################################################
    #Fit the model with overdispersion -- the initial fit
    ######################################################################
    model <- algo.farrington.fitGLM(response,wtime,timeTrend=control$trend,reweight=control$reweight)

    #Stupid check to pass on NULL values from the algo.farrington.fitGLM proc.
    if (is.null(model)) return(model)

    ######################################################################
    #Time trend
    #
    #Check whether to include time trend, to do this we need to check whether
    #1) wtime is signifcant at the 95lvl
    #2) the predicted value is not larger than any observed value
    #3) the historical data span at least 3 years.
    doTrend <- control$trend
    if (control$trend) {
      #is the p-value for the trend significant (0.05) level
      p <- summary.glm(model)$coefficients["wtime",4]
      significant <- (p < 0.05)
      #prediction for time k
      mu0Hat <- predict.glm(model,data.frame(wtime=c(k)),type="response")
      #have to use at least three years of data to allow for a trend
      atLeastThreeYears <- (control$b>=3)
      #no horrible predictions
      noExtrapolation <- mu0Hat <= max(response)
     
      #All 3 criteria have to be met in order to include the trend. Otherwise
      #it is removed. Only necessary to check this if a trend is requested.
      if (!(atLeastThreeYears && significant && noExtrapolation)) {
        doTrend <- FALSE
        model <- algo.farrington.fitGLM(response,wtime,timeTrend=FALSE,reweight=control$reweight)
      }
    } else {
      doTrend <- FALSE
    }
    #done with time trend
    ######################################################################
    
    ######################################################################
    # Calculate prediction & confidence interval                         #
    ######################################################################
    #Predict value - note that the se is the mean CI
    #and not the prediction error of a single observation
    pred <- predict.glm(model,data.frame(wtime=c(k)),dispersion=model$phi,
                        type="response",se.fit=TRUE)
    #Calculate lower and upper threshold
    lu <- algo.farrington.threshold(pred,model$phi,skewness.transform="2/3",alpha=control$alpha)

    ######################################################################
    # If requested show a plot of the fit.
    ######################################################################
    if (control$plot) {
      #Compute all predictions
      data <- data.frame(wtime=c(wtime,k))
      preds <- predict(model,data,,type="response",dispersion=model$phi)

      #Show a plot of the model fit.
      plot(c(wtime, k), c(response,observed[k]),ylim=range(c(observed[data$wtime],lu)),,xlab="time",ylab="No. of counts",main=paste("Prediction at time t=",k))
      #Add the observed value and the upper threshold.
      lines(data$wtime,preds,col=2)
      points(k,lu[2],cex=1,pch=2,col=2)
    }


    ######################################################################
    #Postprocessing steps
    ######################################################################

    #Compute exceedance score unless less than 5 reports during last 4 weeks.
    enoughCases <- (sum(observed[(k-control$limit54[2]):(k-1)])>=control$limit54[1])

    #18 May 2006: Bug/unexpected feature found by Y. Le Strat. 
    #the okHistory variable meant to protect against zero count problems,
    #but instead it resulted in exceedance score == 0 for low counts. 
    #Now removed to be concordant with the Farrington 1996 paper.
    #REMOVEDokHistory <- (pred$fit>1)
    #X <- ifelse(enoughCases && okHistory,
    #            (observed[k] - pred$fit) / (max(lu) - pred$fit),0)
    X <- ifelse(enoughCases,(observed[k] - pred$fit) / (max(lu) - pred$fit),0)

    #Do we have an alarm -- i.e. is observation beyond CI??
    #upperbound only relevant if we can have an alarm (enoughCases)
    trend[k-min(control$range)+1] <- doTrend
    alarm[k-min(control$range)+1] <- (X>1)
    upperbound[k-min(control$range)+1] <- ifelse(enoughCases,lu[2],0)
  }#done looping over all time points

  #Add name and data name to control object.
  control$name <- paste("farrington(",control$w,",",0,",",control$b,")",sep="")
  control$data <- paste(deparse(substitute(disProgObj)))

  # return alarm and upperbound vectors 
  result <- list(alarm = alarm, upperbound = upperbound, trend=trend, 
                 disProgObj=disProgObj, control=control) 
  class(result) <- "survRes" 

  #Done
  return(result) 
}



