################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Scoring rules as discussed in:
### Predictive model assessment for count data
### Czado, C., Gneiting, T. & Held, L. (2009)
### Biometrics 65:1254-1261
###
### Copyright (C) 2010-2012 Michaela Paul
### $Revision: 693 $
### $Date: 2014-01-07 15:21:33 +0100 (Tue, 07 Jan 2014) $
################################################################################


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


## returns logs, rps,ses and dss in reversed!! order
## i.e. the scores for time points n, n-1, n-2,...
scores <- function(object, unit=NULL,sign=FALSE, individual=FALSE)
{
    mu <- object$pred
    x <- object$observed
    size <- object$psi

    if (!is.null(size)) { # NegBin
        size <- exp(size) # transform to parameterization suitable for dnbinom()
        if (ncol(size) != ncol(x)) { # => ncol(size)=1, unit-independent psi
            ## replicate to obtain nrow(size) x nUnits matrix
            size <- matrix(size, nrow=nrow(size), ncol=ncol(x), byrow=FALSE)
        }
    }
    
    if(!is.null(unit)){
        x <- as.matrix(x[,unit])
        mu <- as.matrix(mu[,unit])
        size <- size[,unit]
    }
    
    signXmMu <- if(sign) sign(x-mu) else NULL

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

