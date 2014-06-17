\encoding{latin1}
\name{hhh4_validation}
\alias{oneStepAhead}
\alias{scores}

\title{
 Predictive Model Assessment for \code{hhh4} models
}
\description{
 The function \code{oneStepAhead} computes successive one-step-ahead predictions for
 a (random effects) HHH model fitted by \code{\link{hhh4}}. 

 The function \code{scores} computes a number of (strictly) proper scoring rules based
 on the one-step-ahead predictions.

 See Paul and Held (2011) for further details.
}

\usage{
oneStepAhead(result, tp, type = c("rolling", "first", "final"),
             which.start = c("current", "final"),
             keep.estimates = FALSE, verbose = TRUE, cores = 1)

scores(object, unit = NULL, sign = FALSE, individual = FALSE)
}

\arguments{
  \item{result}{fitted \code{\link{hhh4}} model (class \code{"hhh4"}).}
  \item{tp}{numeric vector of length 1 or 2.
    One-step-ahead predictions are computed from time points 
    \code{tp[1]}, \ldots, \code{tp[2]} (yielding predictions for time
    points \code{tp[1]+1}, ...), where
    \code{tp[2]} defaults to the penultimate time point of
    \code{result$control$subset}.
  }
  \item{type}{The default \code{"rolling"} procedure sequentially 
    refits the model up to each time point in \code{tp} and computes
    the one-step-ahead predictions for the respective next time point.
    The alternative \code{type}s are no true one-step-ahead predictions
    but much faster:
    \code{"first"} will refit the model for the first time point
    \code{tp[1]} only and use this specific fit to calculate all
    subsequent predictions, whereas
    \code{"final"} will just use \code{result} to calculate these.
    The latter case thus gives nothing else than a subset of
    \code{result$fitted.values}, if the \code{tp}'s are part of the
    fitted subset \code{result$control$subset}.}
  \item{which.start}{
    Which initial values should be used when successively
    refitting the model for subsets of the data (up to time point
    \code{tp[1]}, up to \code{tp[1]+1}, ...) if \code{type="rolling"}?
    Default (\code{"current"}) is to use the fitted parameters from the
    previous time point. Alternatively, \code{"final"} means to always
    use the fitted values from \code{result} as initial values for the
    model update. \code{which.start} is ignored for other \code{type}s.
  }
  \item{keep.estimates}{
    logical indicating if parameter estimates and log-likelihoods from
    the successive fits should be returned.
  }
  \item{verbose}{
    non-negative integer (usually in the range \code{0:3}) specifying
    the amount of tracing information to output.
    During \code{hhh4} model updates, the following verbosity is used:
    \code{0} if \code{cores > 1}, otherwise \code{verbose-1} if there
    is more than one time point to predict, otherwise \code{verbose}.
  }
  \item{cores}{the number of cores to use when computing
    the predictions for the set of time points \code{tp} in parallel
    (with \code{\link[parallel]{mclapply}}).
    Note that parallelization is not possible in the default setting
    \code{type="rolling"} and \code{which.start="current"} (use
    \code{which.start="final"} for this to work).}
  \item{object}{result of \code{oneStepAhead}.}
  \item{unit}{integer specifying a specific unit for which the scores
    are computed. If \code{NULL} all units are considered.}
  \item{sign}{logical indicating if the sign of the difference between
    the observed counts and corresponding predictions should be
    returned.}
  \item{individual}{logical indicating if individual scores for all
    units or a mean score over all units should be returned.}
}

\value{
\code{oneStepAhead} returns a list with
\item{pred}{matrix of one-step-ahead predictions}
\item{observed}{matrix with observed counts at the predicted time-points}
\item{psi}{overdispersion parameter on -log-scale in case of a
  negative-binomial model and \code{NULL} otherwise}
\item{allConverged}{logical indicating if all successive fits converged}
If \code{keep.estimates=TRUE}, there are the following additional elements:
\item{coefficients}{matrix with estimated regression parameters from the successive fits}
\item{Sigma.orig}{matrix with estimated variance parameters from the successive fits}
\item{logliks}{matrix with columns \code{"loglikelihood"} and
  \code{"margll"} with their obvious meanings}

\code{scores} returns a matrix with the computed logarithmic, ranked
probability, squared error, Dawid-Sebastiani, and normalized squared
error score as columns.
}

\references{
  Czado, C., Gneiting, T. and Held, L. (2009) Predictive model 
    assessment for count data. Biometrics, \bold{65}, 1254--1261

  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136
}

\author{
Michaela Paul and Sebastian Meyer
}

\seealso{
\code{vignette{"hhh4"}} and \code{\link{hhh4}}.
}

\examples{

## univariate salmonella agona data
data(salmonella.agona)
# convert to sts class
salmonella <- disProg2sts(salmonella.agona)

# generate formula for temporal and seasonal trends
f.end <- addSeason2formula(f = ~ 1 + t, S=1, period=52)
model1 <- list(ar = list(f = ~ 1), end = list(f =f.end),
			  family = "NegBin1")
# run model
result <- hhh4(salmonella, model1)

# do one-step-ahead predictions for the last 5 weeks
pred <- oneStepAhead(result, nrow(salmonella)-5,
                     which.start="final", verbose=FALSE)

\dontshow{
## test equivalence of parallelized version
if (.Platform$OS.type == "unix" && isTRUE(parallel::detectCores() > 1))
    stopifnot(identical(pred,
        oneStepAhead(result, nrow(salmonella)-5,
                     which.start="final", verbose=FALSE, cores=2)))
}

# compute mean scores
colMeans(scores(pred))

# type="final" just means fitted values
stopifnot(identical(
    unname(oneStepAhead(result, nrow(salmonella)-5,
                        type="final", verbose=FALSE)$pred),
    unname(tail(fitted(result), 5))))


\dontrun{
######################################################################
# Do one-step-ahead predictions for the models from the Paul & Held 
# (2011) paper for the influenza data from Bavaria and Baden-Wuerttemberg 
# (this takes some time!)
######################################################################

## see ?hhh4 for a specification of the models

## do 1-step ahead predictions for the last two years

tp <- nrow(fluBYBW)-2*52

val_A0 <- oneStepAhead(res_A0,tp=tp) 
val_B0 <- oneStepAhead(res_B0,tp=tp) 
val_C0 <- oneStepAhead(res_C0,tp=tp) 

val_A1 <- oneStepAhead(res_A1,tp=tp) 
val_B1 <- oneStepAhead(res_B1,tp=tp) 
val_C1 <- oneStepAhead(res_C1,tp=tp) 

val_A2 <- oneStepAhead(res_A2,tp=tp) 
val_B2 <- oneStepAhead(res_B2,tp=tp) 
val_C2 <- oneStepAhead(res_C2,tp=tp) 

val_D <- oneStepAhead(res_D,tp=tp) 


##################################
## compute scores
################################

#scores
vals <- ls(pattern="val_")

nam <- substring(vals,first=5,last=6)

uni <- NULL
indiv <- TRUE

scores_i <- list()
meanScores <- NULL
for(i in seq(along.with=vals)){
  sc <- scores(get(vals[i]),unit=uni, individual=indiv)
  scores_i[[i]] <- sc
  meanScores <- rbind(meanScores,colMeans(sc))
}

names(scores_i) <- nam
rownames(meanScores) <- nam

##comparison with best model B2 

compareWithBest <- function(best, whichModels, whichScores=1:3, nPermut=9999, seed=1234){
  set.seed(seed)
  pVals <- NULL
  for(score in whichScores){
    p <- c()
    for(model in whichModels){
      if(model==best) p <- c(p,NA)
      else p <- c(p,permutationTest(scores_i[[model]][,score],scores_i[[best]][,score],
     plot=TRUE,nPermutation=nPermut)$pVal.permut)
    }  
    pVals <- cbind(pVals,p)
  }
  return(pVals)
}


pVals_flu <- compareWithBest(best=6, whichModels=1:10,seed=2059710987)
rownames(pVals_flu) <- nam

}

}
\keyword{ts}
\keyword{regression}
