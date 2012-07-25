#####################################################
## Scoring rules as discussed in:
## Predictive model assessment for count data
## Czado, C., Gneiting, T. & Held, L. (2009)
## Biometrics 65:1254-1261
######################################################

## logarithmic score
# logs(P,x) = - log(P(X=x))
logScore <- function(x,mu, size=NULL){
	if(is.null(size))
		- dpois(x,lambda=mu,log=TRUE)
	else - dnbinom(x, mu=mu,size=size,log=TRUE)
}

## squared error score
# ses(P,x) =(x-mu_p)^2
ses <- function(x, mu){
  (x-mu)^2
}

## normalized squared error score
# nses(P,x) =((x-mu_p)/sigma_p)^2
nses <- function(x, mu,size=NULL){
  if(!is.null(size)){
    sigma2 <- mu*(1+mu/size)
  } else sigma2 <- mu
  
  ((x-mu)^2)/sigma2
}

## Dawid-Sebastiani score
# dss(P,x) = ((x-mu_p)/sigma_p)^2 + 2*log(sigma_p)
dss <- function(x,mu,size=NULL){
  if(!is.null(size)){
    sigma2 <- mu*(1+mu/size)
  } else sigma2 <- mu
  
  ((x-mu)^2)/sigma2 +log(sigma2)
}

## ranked probability score
# rps(P,x) =sum_0^Kmax { P(X<=k) - 1(x <=k)}^2
rps.one <- function(x, mu,size=NULL,k=40,eps=1e-10){
	# determine variance of distribution
	if(is.null(size)){
	  se <- sqrt(mu)
	} else se <- sqrt(mu*(1+mu/size))
	
	# determine the maximum number of summands as Kmax= mean+k*se
	kmax <- ceiling(mu + k*se)
	
	# compute 1(x <=k)
	ind <- 1*(x < (1:(kmax+1)))
	
	# compute P(X<=k)
	# Poisson case
	if(is.null(size)){
		Px <- ppois(0:kmax,lambda=mu)
	} else  Px <- pnbinom(0:kmax,mu=mu,size=size)
	
	#determine precision
	if((1-tail(Px,1))^2 > eps)
	 cat("precision of finite sum not smaller than ", eps,"\n")
	
	
	# compute rps
	sum((Px-ind)^2)
}

rps <- function(x,mu,size=NULL,k=40){
	n <- length(x)
	if(length(mu)==1)
		mu <- rep(mu,n)
	if(!is.null(size) & length(size)==1)
		size <- rep(size,n)
	
	res <- sapply(1:n, function(i) rps.one(x=x[i],mu=mu[i],size=size[i],k=k) )
	matrix(res,ncol=ncol(as.matrix(x)),byrow=FALSE)
}

##
## returns logs, rps,ses and dss in reversed!! order
## i.e. the scores for time points n, n-1, n-2,...
scores <- function(object, unit=NULL,sign=FALSE, individual=FALSE){
	mu <- object$mean
	size <- object$psi
	x <- object$x
	
	if(!is.null(unit)){
		x <- as.matrix(x[,unit])
		mu <- as.matrix(mu[,unit])
		if(ncol(size)>1){
			size <- as.matrix(size[,unit])
		} 
	}
	
	#negbin or poisson?
	if(any(is.na(size))){
		size <- NULL
	} else if(ncol(size)!=ncol(x)){
		size <- matrix(exp(size), nrow=nrow(size), ncol=ncol(x),byrow=FALSE)
	} else {
      size <- exp(size)  
   }
	
	if(sign)
		signXmMu <- sign(x-mu)
	else signXmMu <- NULL

	#compute average scores for unit
	log.score <- apply(as.matrix(logScore(x=x,mu=mu,size=size)),MARGIN=2,rev)
	rp.score <- apply(as.matrix(rps(x=x,mu=mu,size=size)),MARGIN=2,rev)
	se.score <- apply(as.matrix(ses(x=x,mu=mu)), MARGIN=2, rev)
	nse.score <- apply(as.matrix(nses(x=x,mu=mu,size=size)),MARGIN=2,rev)
	ds.score <- apply(as.matrix(dss(x=x, mu=mu, size=size)), MARGIN=2,rev)
	
	if(is.null(unit)){
      if(individual){
        log.score <- c(log.score)
        rp.score <- c(rp.score)
        se.score <- c(se.score)
        nse.score <- c(nse.score)
        ds.score <- c(ds.score)
      } else {
        log.score <- rowMeans(log.score)
        rp.score <- rowMeans(rp.score)
        se.score <- rowMeans(se.score)
        nse.score <- rowMeans(nse.score)
        ds.score <- rowMeans(ds.score)
      }    
	}
	
	result <- cbind(logs=log.score,rps=rp.score,ses=se.score,dss=ds.score,nses=nse.score,signXmMu=signXmMu)
	return(result)
}


#####################################################
## non-randomized version of the PIT histogram
##
## Params: 
## J - number of bins 
## x - observed data
## pdistr - cumulative distribution function, either pnbinom or ppois
## ... - arguments for pdistr (i.e. lambda for ppois, mu and size for pnbinom)
####################################################
pit.one <- function(u,x,pdistr=pnbinom,...){

   Px <- pdistr(x,...)
   Pxm1 <- pdistr(x-1,...)
   F_u <- ifelse(u< Pxm1 | u== Pxm1, 0, pmin(1,(u-Pxm1)/(Px-Pxm1) ) )
   F_u[x==0] <- pmin(1,u/pdistr(0,...))[x==0]
   if(u==1){
    F_u <- 1
   } else if(u==0){
    F_u <- 0  
   }

   return(mean(F_u))
}


pit <- function(J=10,x,pdistr=pnbinom,...){
  F_u.bar <- sapply((0:J)/J,pit.one,x=x,pdistr=pdistr,...)
  f_j <- J*diff(F_u.bar)
  
  erg <- list(breaks=(0:J)/J,counts=f_j, density=f_j,mids=(0:(J-1))/J+diff((0:J)/J)/2,
             xname="PIT",equidist=TRUE)
  class(erg) <- "histogram"
  return(erg)
}

#######################################################
## scores1, scores2 - vector with scores from two models
#########################################################
permutationTest <- function(score1,score2, nPermutation=9999,plot=FALSE){
  meanScore1 <- mean(score1)
  meanScore2 <- mean(score2)
  diffObserved <- meanScore1 -meanScore2
  
  nTime <- length(score1)
  diffMean <- rep(NA,nPermutation)
  
  for(i in 1:nPermutation){
    sel <- rbinom(nTime, size=1, prob=0.5)
    g1 <- (sum(score1[sel==0]) + sum(score2[sel==1]))/nTime
    g2 <- (sum(score1[sel==1]) + sum(score2[sel==0]))/nTime
    diffMean[i] <- g1-g2
  }

  if(plot){
    hist(diffMean, nclass=50, prob=TRUE,xlab="Difference between means",main="")
    abline(v=diffObserved,col=4)
  }
  
  pVal <- (1+sum(abs(diffMean)>=abs(diffObserved)))/(nPermutation+1)

  pTtest <- t.test(score1,score2,paired=TRUE)$p.value
  
  cat("mean difference=",diffObserved,"\tp(permutation) =",pVal,"\tp(paired t-test) =",pTtest,"\n")
  return(list(diffObs=diffObserved, pVal.permut=pVal,pVal.t=pTtest))
}

