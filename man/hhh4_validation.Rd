\encoding{latin1}
\name{hhh4_validation}
\alias{oneStepAhead}
\alias{scores}

\title{
 Predictive model assessment for a HHH4 model
}
\description{
 The function \code{oneStepAhead} computes successive one-step-ahead predictions for
 a (random effects) HHH model resulting from a call to \code{\link{hhh4}}. 

 The function \code{scores} computes a number of (strictly) proper scoring rules based
 on the one-step-ahead predictions.

 See Paul and Held (2011) for further details.
}
\usage{
oneStepAhead(result, tp, verbose = TRUE, keep.estimates = FALSE)

scores(object, unit = NULL, sign = FALSE, individual = FALSE)
}
\arguments{
  \item{result}{Result of a call to \code{\link{hhh4}}.
}
  \item{tp}{
  One-step-ahead predictions are computed for time-points 
  \code{tp+1}, \ldots, \code{nrow(result$stsObj)}.
}
  \item{verbose}{
   Logical if additional information is to be printed during the computations.
}
  \item{keep.estimates}{
  Logical if parameter estimates from the successive fits should be returned.
}
  \item{object}{Result of a call to \code{oneStepAhead}.
}
  \item{unit}{Integer specifying a specific unit for which the scores are computed. If \code{null} all units are considered.
}
  \item{sign}{Logical indication if the sign of the difference between the observed counts and corresponding predictions should be returned.
}
  \item{individual}{Logical indicating if individual scores for all units or a mean score over all units should be returned.
}
}
%\details{
%}
\value{
\code{oneStepAhead} returns a list with
 \item{mean}{matrix with one-step-ahead predictions}
 \item{psi}{overdispersion parameter on log scale in case of a negative binomial model and \code{NULL} otherwise}
 \item{x}{matrix with observed counts at time-points 
  \code{tp+1}, \ldots, \code{nrow(result$stsObj)} }
 \item{allConverged}{Logical indicating if all successive fits converged}
 \item{params}{matrix with estimated regression parameters from the successive fits if
  \code{keep.estimates = TRUE}, and \code{NULL} otherwise.}
 \item{variances}{matrix with estimated variance components from the successive fits if
  the model contains random effects and \code{keep.estimates = TRUE},  and \code{NULL} otherwise.}

\code{scores} returns a matrix with the computed logarithmic, ranked probability, squared error, Dawid-Sebastiani, and normalized squared error score as columns.
}
\references{
  Czado, C., Gneiting, T. and Held, L. (2009) Predictive model 
    assessment for count data. Biometrics, \bold{65}, 1254--1261

  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136
}
\author{
M. Paul
}
\seealso{
See the vignette \code{vignette{"hhh4"}} and \code{\link{hhh4}}.
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
