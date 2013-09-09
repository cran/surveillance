\encoding{latin1}
\name{hhh4_validation}
\alias{oneStepAhead}
\alias{scores}

\title{
 Predictive Model Assessment for HHH4 models
}
\description{
 The function \code{oneStepAhead} computes successive one-step-ahead predictions for
 a (random effects) HHH model fitted by \code{\link{hhh4}}. 

 The function \code{scores} computes a number of (strictly) proper scoring rules based
 on the one-step-ahead predictions.

 See Paul and Held (2011) for further details.
}

\usage{
oneStepAhead(result, tp, which.start = c("current", "final", "none"),
             keep.estimates = FALSE, verbose = TRUE)

scores(object, unit = NULL, sign = FALSE, individual = FALSE)
}

\arguments{
  \item{result}{Fitted \code{\link{hhh4}} model object of class \code{"ah4"}.}
  \item{tp}{numeric vector of length 1 or 2.
    One-step-ahead predictions are computed from time points 
    \code{tp[1]}, \ldots, \code{tp[2]} (yielding predictions for time
    points \code{tp[1]+1}, ...), where
    \code{tp[2]} defaults to \code{nrow(result$stsObj)}.
  }
  \item{which.start}{
    Which initial values should be used when successively
    refitting the model for subsets of the data (up to time point
    \code{tp[1]}, up to \code{tp[1]+1}, ...)?
    Default (\code{"current"}) is to use the fitted parameters from the
    previous time point. Alternatively, \code{"final"} means to always
    use the fitted values from \code{result} as initial values for the
    model update, and \code{"none"} means to use \code{result} all along
    time (no updating at all). The latter case runs much faster but reduces
    to comparing the fitted means from \code{result} with the observed
    data (no one-step-ahead predictions).
  }
  \item{keep.estimates}{
    logical indicating if parameter estimates and log-likelihoods from
    the successive fits should be returned.
  }
  \item{verbose}{
    logical indicating if additional information is to be printed during
    the computations. For the sequential \code{hhh4} model updates, a
    verbosity of \code{verbose-1} is used if there is more than one time
    point to predict, otherwise \code{verbose}.
  }
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
\item{pred}{matrix with one-step-ahead predictions}
\item{psi}{overdispersion parameter on log-scale in case of a negative
  binomial model and \code{NULL} otherwise}
\item{observed}{matrix with observed counts at the predicted time-points}
\item{allConverged}{logical indicating if all successive fits converged}
If \code{keep.estimates = TRUE}, the following elements will be non-\code{NULL}:
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
pred <- oneStepAhead(result, nrow(salmonella)-5)

# compute mean scores
colMeans(scores(pred))

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
