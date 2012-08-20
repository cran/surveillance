################################################################################
# Authors: Sebastian Meyer, Michael Hoehle
# Date: 4 Jun 2009
#
# This file contains methods for generic functions for objects of class
# "twinSIR", specifically:
# - coef and vcov: enabling the use of function confint to calculate Wald
#                  confidence intervals for the parameter estimates.
# - logLik: enables the use of function AIC
# - AIC, extractAIC: compute AIC or OSAIC depending on argument 'one.sided'
# - print, summary, print.summary, plot (intensityPlot), ...
# - profile: Calculates the profile log-likelihood (-> likelihood.ci)
################################################################################

### don't need a specific coef-method (identical to stats:::coef.default)
## coef.twinSIR <- function (object, ...)
## {
##     object$coefficients
## }

# asymptotic variance-covariance matrix (inverse of fisher information matrix)
vcov.twinSIR <- function (object, ...)
{
    solve(object$fisherinfo)
}

logLik.twinSIR <- function (object, ...)
{
    r <- object$loglik
    attr(r, "df") <- length(coef(object))
    class(r) <- "logLik"
    r
}

# Note: pz is determined by scanning the names of coef(object),
#       thus the 'model' component is not necessary
# See the Hughes and King (2003) paper for details
.OSAICpenalty <- function (twinSIRobject, k = 2, nsim = 1e3)
{
    theta <- coef(twinSIRobject)
    npar <- length(theta)
    pz <- length(grep("cox\\([^)]+\\)", names(theta), ignore.case = FALSE,
                      perl = FALSE, fixed = FALSE, useBytes = FALSE,
                      invert = FALSE))
    px <- npar - pz   # number of constrained (non-negative) parameters
    
    penalty <- if (px == 0L) {
        k * pz   # default AIC penalty (with k = 2)
    } else if (px == 1L) {
        k * (pz + 0.5)
    } else if (px == 2L) {
        Sigma <- vcov(twinSIRobject)   # parameter covariance matrix
        rho <- cov2cor(Sigma[1:2,1:2])[1,2]
        as <- acos(rho)/2/pi
        w <- c(as, 0.5, 0.5-as)
        k * sum(w * (pz + 0:2))   # = k * sum(w * (npar - px + 0:2))
    } else { # px > 2
        cat("Computing OSAIC weights for px =",px,"epidemic covariates",
            "based on", nsim, "simulations...\n")
        W <- vcov(twinSIRobject)[1:px,1:px]
        #The simulation approach has unresolved problems for some situations (initial constraints
        #can fail). Catch this operation until exact cause has been investigated
        tryCatch( { 
          w.sim <- w.chibarsq.sim(p=px, W=W, N=nsim)
          #c.f. (12) in Hughes & King (2003), r_i=px, m=0:px, ki=npar
          #as npar=pz+px, we have that npar-px = pz, hence the sum is
          k * sum(w.sim * (pz + 0:px))
        }, error=function(e) NA)
    }
    
    attr(penalty, "exact") <- px <= 2
    penalty
}

AIC.twinSIR <- function (object, ..., k = 2, one.sided = NULL, nsim = 1e3)
{
    AIC.default <- match.call()
    AIC.default$one.sided <- NULL
    AIC.default$nsim <- NULL
    AIC.default[[1]] <- call(":::", as.name("stats"), as.name("AIC.default"))
    
    if (is.null(one.sided)) {
        one.sided <- object$method == "L-BFGS-B"
    }
    
    res <- if (!one.sided) {
        eval(AIC.default, parent.frame())
    } else {
        penalty <- .OSAICpenalty(object, k = k, nsim = nsim)
        edf <- length(coef(object))
        AIC.default$k <- penalty/edf
        eval(AIC.default, parent.frame())
    }
    
    attr(res, "type") <- if (one.sided) "One-sided AIC" else "Standard AIC"
    attr(res, "exact") <- if (one.sided) attr(penalty, "exact") else TRUE
    res
}

extractAIC.twinSIR <- function (fit, scale = 0, k = 2, one.sided = NULL,
    nsim = 1e3, ...)
{
    if (is.null(one.sided)) {
        one.sided <- fit$method == "L-BFGS-B"
    }
    
    loglik <- logLik(fit)
    edf <- attr(loglik, "df")
    penalty <- if (one.sided) {
                   .OSAICpenalty(fit, k = k, nsim = nsim)   # one-sided AIC
               } else {
                   k * edf                                  # default AIC
               }
    res <- c(edf = edf, AIC = -2 * c(loglik) + penalty)            
    
    attr(res, "type") <- if (one.sided) "One-sided AIC" else "Standard AIC"
    attr(res, "exact") <- if (one.sided) attr(penalty, "exact") else TRUE
    res
}

print.twinSIR <- function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    print.default(format(coef(x), digits=digits), print.gap = 2, quote = FALSE)
    cat("\nLog-likelihood: ", format(logLik(x), digits=digits), "\n", sep = "")
    if (!x$converged) {
        cat("\nWARNING: OPTIMIZATION DID NOT CONVERGE!\n")
    }
    cat("\n")
    invisible(x)
}

summary.twinSIR <- function (object,
    correlation = FALSE, symbolic.cor = FALSE, ...)
{
    ans <- object[c("call", "converged", "counts", "intervals", "nEvents")]
    ans$cov <- vcov(object)
    est <- coef(object)
    se <- sqrt(diag(ans$cov))
    zval <- est/se
    pval <- 2 * pnorm(abs(zval), lower.tail = FALSE)
    ans$coefficients <- cbind(est, se, zval, pval)
    dimnames(ans$coefficients) <- list(names(est),
        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    if (correlation) {
        ans$correlation <- cov2cor(ans$cov)
        ans$symbolic.cor <- symbolic.cor
    }
    ans$loglik <- logLik(object)
    aic <- extractAIC(object, ...)
    ans$aic <- as.vector(aic[2L])   # remove 'edf' element
    attributes(ans$aic) <- attributes(aic)[c("type", "exact")]
    class(ans) <- "summary.twinSIR"
    ans
}

print.summary.twinSIR <- function (x,
    digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
    signif.stars = getOption("show.signif.stars"), ...)
{
    cat("\nCall:\n")
    print.default(x$call)
    cat("\nCoefficients:\n")
    coefs <- x$coefficients
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars,
                 na.print = "NA", ...)
    nEvents <- x$nEvents
    nh0 <- length(nEvents)
    if (nh0 < 2L) {
        cat("\nTotal number of infections: ", nEvents, "\n")
    } else {
        cat("\nBaseline intervals:\n")
        intervals <- character(nh0)
        for(i in seq_len(nh0)) {
            intervals[i] <-
            paste("(",
                  paste(format(x$intervals[c(i,i+1L)],trim=TRUE), collapse=";"),
                  "]", sep = "")
        }
        names(intervals) <- paste("logbaseline", seq_len(nh0), sep=".")
        print.default(rbind("Time interval" = intervals,
                            "Number of events" = nEvents),
                      quote = FALSE, print.gap = 2)
    }
    cat("\n", attr(x$aic, "type"), ": ", format(x$aic, digits=max(4, digits+1)),
        if (!attr(x$aic, "exact")) "\t(simulated penalty weights)" else "",
        sep = "")
    cat("\nLog-likelihood:", format(x$loglik, digits = digits))
    cat("\nNumber of log-likelihood evaluations:", x$counts[1], "\n")
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1L) {
        cat("\nCorrelation of Coefficients:\n")
        if (is.logical(symbolic.cor) && symbolic.cor) {
            correl <- symnum(correl, abbr.colnames = NULL)
            correlcodes <- attr(correl, "legend")
            attr(correl, "legend") <- NULL
            print(correl)
            cat("---\nCorr. codes:  ", correlcodes, "\n", sep="")
        } else {
            correl <- format(round(correl, 2), nsmall = 2, digits = digits)
            correl[!lower.tri(correl)] <- ""
            print(correl[-1, -p, drop = FALSE], quote = FALSE)
        }
        }
    }
    if (!x$converged) {
        cat("\nWARNING: OPTIMIZATION DID NOT CONVERGE!\n")
    }
    cat("\n")
    invisible(x)
}


### Plot method for twinSIR (wrapper for intensityplot)

plot.twinSIR <- function (x, which, ...) # defaults for 'which' are set below
{
    cl <- match.call()
    cl[[1]] <- as.name("intensityplot")
    eval(cl, envir = parent.frame())
}

formals(plot.twinSIR)$which <- formals(intensityplot.twinSIR)$which




######################################################################
# Function to compute likelihood based confidence interval, basically
# the two solutions to
#            f(\theta) = l(\theta)-l(\hat{theta)) + 1/2 dchisq(1-alpha,df=1)=0
# are found.
#
#
# Parameters:
#  logliktilde - normalized likelihood function(theta, ...)
#  theta.hat - the MLE
#  lower - search interval [lower,theta.hat] for f=0
#  upper - search interval [theta.hat,upper] for f=0
#  alpha - confidence level (see Equation 2.6 in Pawitan (2003)
#  ... - additional arguments passed to function logliktilde
######################################################################

likelihood.ci <- function (logliktilde, theta.hat, lower, upper,
    alpha = 0.05, ...)
{
  # Highest Likelihood intervall -- target function
  f <- function(theta, ...) { 
    logliktilde(theta, ...) + 1/2*qchisq(1-alpha, df=1)
  }
  # Compute upper and lower boundary numerically
  hl.lower <- uniroot(f, interval = c(lower, theta.hat), ...)$root
  hl.upper <- uniroot(f, interval = c(theta.hat, upper), ...)$root
  return(c(hl.lower,hl.upper))
}


######################################################################
# Function to compute estimated and profile likelihood based
# confidence intervals. Heavy computations might be necessary!
#
#Params:
# fitted - output from a fit with twinSIR
# profile - list with 4D vector as entries - format:
#               c(index, lower, upper, grid size)
#           where index is the index in the coef vector
#                 lower and upper are the parameter limits (can be NA)
#                 grid size is the grid size of the equally spaced grid
#                 between lower and upper (can be 0)
# alpha - (1-alpha)% profile likelihood CIs are computed.
#         If alpha <= 0 then no CIs are computed
# control - control object to use for optim in the profile loglik computations
#
# Returns:
#  list with profile loglikelihood evaluations on the grid
#  and highest likelihood and wald confidence intervals
######################################################################

profile.twinSIR <- function (fitted, profile, alpha = 0.05,
    control = list(fnscale = -1, factr = 1e1, maxit = 100), ...)
{
  ## Check that input is ok
  profile <- as.list(profile)
  if (length(profile) == 0L) {
    stop("nothing to do")
  }
  lapply(profile, function(one) {
    if (length(one) != 4L) {
      stop("each profile entry has to be of form ",
           "'c(index, lower, upper, grid size)'")
    }})
  if (is.null(fitted[["model"]])) {
    stop("'fitted' must contain the model component")
  }
  
  px <- ncol(fitted$model$X)
  pz <- ncol(fitted$model$Z)
  
  ## Control of the optim procedure
  if (is.null(control[["fnscale",exact=TRUE]])) { control$fnscale <- -1 }
  if (is.null(control[["factr",exact=TRUE]])) { control$factr <- 1e1 }
  if (is.null(control[["maxit",exact=TRUE]])) { control$maxit <- 100 }

  
  ## Estimated normalized likelihood function
  ltildeestim <- function(thetai,i) {
    theta <- theta.ml
    theta[i] <- thetai
    with(fitted$model,
      .loglik(theta, X=X, Z=Z, survs=survs, weights=weights)) - loglik.theta.ml
  }

  ## Profile normalized likelihood function
  ltildeprofile <- function(thetai,i)
  {
    emptyTheta <- rep(0, length(theta.ml))
      
    # Likelihood l(theta_{-i}) = l(theta_i, theta_i)
    ltildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      with(fitted$model,
       .loglik(theta, X=X, Z=Z, survs=survs, weights=weights)) - loglik.theta.ml
    }
    # Score function of all params except thetaminusi
    stildethetaminusi <- function(thetaminusi) {
      theta <- emptyTheta
      theta[-i] <- thetaminusi
      theta[i] <- thetai
      with(fitted$model,
        .score(theta, X=X, Z=Z, survs=survs, weights=weights))[-i]
    }
      
    # Call optim using L-BFGS-B. For harder constrains we need constr.Optim
    lower <- if (fitted$method == "L-BFGS-B") {
               c(rep(0,px),rep(-Inf,pz))[-i]
             } else {
               -Inf
             }
    upper <- if (fitted$method == "L-BFGS-B") {
               c(rep(Inf,px),rep(Inf,pz))[-i]
             } else {
               Inf
             }
    resOthers <- tryCatch(with(fitted$model,
            optim(theta.ml[-i], fn = ltildethetaminusi, gr = stildethetaminusi,
                  method = fitted$method, control = control,
                  lower = lower, upper = upper)),
            warning = function(w) print(w), error = function(e) list(value=NA))
    resOthers$value
  }
  
  ## Initialize
  theta.ml <- coef(fitted)
  loglik.theta.ml <- c(logLik(fitted))
  se <- sqrt(diag(vcov(fitted)))
  resProfile <- list()
  
  ## Perform profile computations for all requested parameters
  cat("Evaluating the profile logliks on a grid...\n")
  for (i in 1:length(profile))
  {
    cat("i= ",i,"/",length(profile),"\n")
    #Index of the parameter in the theta vector
    idx <- profile[[i]][1]
    #If no borders are given use those from wald intervals (unconstrained)
    if (is.na(profile[[i]][2])) profile[[i]][2] <- theta.ml[idx] - 3*se[idx]
    if (is.na(profile[[i]][3])) profile[[i]][3] <- theta.ml[idx] + 3*se[idx]
    #Evaluate profile loglik on a grid (if requested)
    if (profile[[i]][4] > 0) {
      thetai.grid <- seq(profile[[i]][2],profile[[i]][3],length=profile[[i]][4])
      resProfile[[i]] <- matrix(NA, nrow = length(thetai.grid), ncol = 4L,
        dimnames = list(NULL, c("grid","profile","estimated","wald")))
      
      for (j in 1:length(thetai.grid)) {
        cat("\tj= ",j,"/",length(thetai.grid),"\n")
        resProfile[[i]][j,] <- c(thetai.grid[j],
           ltildeprofile(thetai.grid[j],idx),
           ltildeestim(thetai.grid[j],idx),
#9 June 2009: Bug discovered by L. Held. as part of paper revision. C.f. Pawitan p.63
           - 1/2*(1/se[idx]^2)*(thetai.grid[j] - theta.ml[idx])^2)
      }
    }
  }
#9 June 2009. This did not work.
#  names(resProfile) <- names(theta.ml)[sapply(profile, function(x) x[4L]) > 0]
   names(resProfile) <- names(theta.ml)[sapply(profile, function(x) x[1L])]
  
  ## Profile likelihood intervals
  lower <- if (fitted$method == "L-BFGS-B") {
             c(rep(0,px),rep(-Inf,pz))
           } else {
             -Inf
           }
  ciProfile <- matrix(NA, nrow = length(profile), ncol = 5L,
    dimnames = list(NULL, c("idx","hl.low","hl.up","wald.low","wald.up")))
  if (alpha > 0) {
    cat("Computing likelihood ratio intervals...\n")
    for (i in 1:length(profile))
    {
      cat(i,"/", length(profile),"\n") 
      #Index of the parameter in the theta vector
      idx <- profile[[i]][1]
      #Compute highest likelihood intervals
      ci.hl <- tryCatch(
        likelihood.ci(ltildeprofile, theta.hat = theta.ml[idx],
                      lower = max(lower[idx], theta.ml[idx]-5*se[idx]),
                      upper = theta.ml[idx]+5*se[idx], alpha = alpha, i = idx),
        warning = function(w) print(w),
        error = function(e) rep(NA,2))
      #Wald intervals based on expected fisher information
      ci.wald <- theta.ml[idx] + c(-1,1) * qnorm(1-alpha/2) * se[idx]
      ciProfile[i,] <- c(idx, ci.hl, ci.wald)
    }
    rownames(ciProfile) <- names(theta.ml)[ciProfile[,1]]
  }
  
  return(list(lp=resProfile, ci.hl=ciProfile, profileObj=profile))
}


######################################################################
# Extract the "residual process" (cf. Ogata, 1988), i.e. the
# fitted cumulative intensity at the event times.
# -> "generalized residuals similar to those discussed in Cox and Snell (1968)"
######################################################################

residuals.twinSIR <- function(object, ...)
{
  #Extract event and stop-times
  eventTimes <- attr(object$model$survs,"eventTimes")
  sortedStop <- sort(unique(object$model$survs[,"stop"]))
  eventTimesIdx <- match(eventTimes, sortedStop)
  
  #Dimensions and zero vector (in case we need it)
  nTimes <- nrow(object$model$X)
  zerovec <- numeric(nTimes)

  # Extract the fitted model params
  px <- ncol(object$model$X)
  pz <- ncol(object$model$Z)
  theta <- coef(object)
  alpha <- theta[seq_len(px)]
  beta <- theta[px+seq_len(pz)]

  # Initialize e, h and thus lambda
  if (px > 0) { e <- as.vector(object$model$X %*% as.matrix(alpha)) } else { e <- zerovec }
  if (pz > 0) { h <- as.vector(exp(object$model$Z %*% as.matrix(beta))) } else { h <- zerovec }
  lambda <- (e + h)

  #Determine bloks
  BLOCK <- as.numeric(factor(object$model$survs$start))

  # lambda_i integrals, i.e. integral of \lambda_i until t for each individual
  dt <- object$model$survs[,"stop"] - object$model$survs[,"start"]

  #Easier - no individual summations as they are all summed anyhow afterwards
  intlambda <- tapply(object$model$weights * lambda* dt, BLOCK, sum)

  #Compute cumulative intensities (Ogata (1988): "residual process")
  tau <- cumsum(intlambda)[eventTimesIdx]
  tau
}

