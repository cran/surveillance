\encoding{latin1}
\name{hhh4_validation}
\alias{oneStepAhead}
\alias{scores}
\alias{scores.oneStepAhead}
\alias{scores.hhh4}
\alias{calibrationTest.oneStepAhead}
\alias{pit.oneStepAhead}
\alias{pit.hhh4}

\title{Predictive Model Assessment for \code{hhh4} Models}
\description{
  The function \code{oneStepAhead} computes successive one-step-ahead
  predictions for a (random effects) HHH model fitted by \code{\link{hhh4}}. 

  The function \code{scores} computes a number of (strictly) proper
  scoring rules based on such one-step-ahead predictions
  (or, for the \code{hhh4}-method, directly based on the fitted values).

  See Paul and Held (2011) for further details.

  There are also \code{\link{calibrationTest}} and \code{\link{pit}}
  methods for \code{oneStepAhead} predictions.
}

\usage{
oneStepAhead(result, tp, type = c("rolling", "first", "final"),
             which.start = c("current", "final"),
             keep.estimates = FALSE, verbose = TRUE, cores = 1)

scores(x, ...)
\method{scores}{oneStepAhead}(x, which = c("logs", "rps", "dss", "ses"),
       units = NULL, sign = FALSE, individual = FALSE, reverse = TRUE, ...)
\method{scores}{hhh4}(x, which = c("logs", "rps", "dss", "ses"),
       subset = x$control$subset, units = seq_len(x$nUnit), sign = FALSE, ...)

\method{calibrationTest}{oneStepAhead}(x, ...)

\method{pit}{oneStepAhead}(x, ...)
\method{pit}{hhh4}(x, subset = x$control$subset, units = seq_len(x$nUnit), ...)
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
    Which initial parameter values should be used when successively
    refitting the model to subsets of the data (up to time point
    \code{tp[1]}, up to \code{tp[1]+1}, ...) if \code{type="rolling"}?
    Default (\code{"current"}) is to use the parameter estimates from the
    previous time point, and \code{"final"} means to always
    use the estimates from \code{result} as initial values.
    Alternatively, \code{which.start} can be a list of \code{start}
    values as expected by \code{\link{hhh4}}, which then replace
    the corresponding estimates from \code{result} as initial values.
    This argument is ignored for \dQuote{non-rolling} \code{type}s.
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
  
  \item{x}{an object of class \code{"oneStepAhead"} or \code{"hhh4"}.}
  \item{which}{character vector determining which scores to compute.
    The package \pkg{surveillance} implements the following proper
    scoring rules: logarithmic score (\code{"logs"}), ranked probability
    score (\code{"rps"}), Dawid-Sebastiani score (\code{"dss"}), and
    squared error score (\code{"ses"}). The normalized SES
    (\code{"nses"}) is also available but it is improper and not
    computed by default.\cr
    It is possible to name own scoring rules in \code{which}. These
    must be functions of \code{(x, mu, size)}, vectorized in all arguments
    (time x unit matrices) except that \code{size} is \code{NULL}
    in case of a Poisson model.
    See the supplied scoring rules for guidance, e.g.,
    \code{\link{logs}} or \code{\link{dss}}.}
  \item{subset}{
    subset of time points for which to calculate the scores
    (or the PIT histogram).
    Defaults to the subset used for fitting the model.}
  \item{units}{integer or character vector indexing the units for which
    the scores (or the PIT histogram) should be computed.
    By default (\code{NULL}) all units are considered.}
  \item{sign}{logical indicating if the function should also return
    \code{sign(x-mu)}, i.e., the sign of the difference between
    the observed counts and corresponding predictions.
    This does not really make sense when averaging over multiple
    \code{units} with \code{individual=FALSE}.}
  \item{individual}{logical indicating if the individual scores of the
    \code{units} should be returned. By default (\code{FALSE}), the
    individual scores are averaged over all \code{units}.}
  \item{reverse}{logical indicating if the rows (time points) should be
    reversed in the result (historical default for the
    \code{oneStepAhead}-method).}
  \item{\dots}{Unused by the \code{scores}-methods.\cr
    For the \code{calibrationTest}-method, further arguments are passed
    to \code{\link{calibrationTest.default}}, e.g., \code{which} to
    select a scoring rule.\cr
    For the \code{pit}-methods, further arguments are passed to
    \code{\link{pit.default}}.}
}

\value{
\code{oneStepAhead} returns a list (of class \code{"oneStepAhead"})
with the following components:
\item{pred}{one-step-ahead predictions in a matrix, where each row
  corresponds to one of the time points requested via the argument
  \code{tp}, and which has \code{ncol(result$stsObj)}
  unit-specific columns. The rownames indicate the predicted time points
  and the column names are identical to \code{colnames(result$stsObj)}.}
\item{observed}{matrix with observed counts at the predicted time
  points. It has the same dimensions and names as \code{pred}.}
\item{psi}{in case of a negative-binomial model, a matrix of the
  estimated overdispersion parameter(s) at each time point on 
  the internal -log-scale (1 column if \code{"NegBin1"},
  \code{ncol(observed)} columns if \code{"NegBinM"} or shared overdispersion). 
  For a \code{"Poisson"} model, this component is \code{NULL}.}
\item{allConverged}{logical indicating if all successive fits
  converged.}

If \code{keep.estimates=TRUE}, there are the following additional elements:
\item{coefficients}{matrix of estimated regression parameters from the successive fits.}
\item{Sigma.orig}{matrix of estimated variance parameters from the successive fits.}
\item{logliks}{matrix with columns \code{"loglikelihood"} and
  \code{"margll"} with their obvious meanings.}


The function \code{scores} computes the scoring rules specified in the
argument \code{which}.
If multiple \code{units} are selected and \code{individual=TRUE}, the
result is an array of dimensions
\code{c(nrow(pred),length(units),5+sign)} (up to \pkg{surveillance}
1.8-0, the first two dimensions were collapsed to give a matrix).
Otherwise, the result is a matrix with \code{nrow(pred)} rows and
\code{5+sign} columns. If there is only one predicted time point, the
first dimension is dropped in both cases.
CAVE: For historical reasons, \code{scores.oneStepAhead} by default
\code{reverse}s the order of the rows (time points)!
The \code{hhh4}-method of \code{scores} does not reverse rows.


The \code{\link{calibrationTest}}- and \code{\link{pit}}-methods are
just convenient wrappers around the respective default methods.
See also \code{\link{calibrationTest.hhh4}}.
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
\code{vignette("hhh4")} and \code{\link{hhh4}}.
}

\examples{
### univariate salmonella agona data

data("salmonella.agona")
## convert to sts class
salmonella <- disProg2sts(salmonella.agona)

## generate formula for temporal and seasonal trends
f.end <- addSeason2formula(f = ~1 + t, S=1, period=52)
model1 <- list(ar = list(f=~1), end = list(f=f.end), family = "NegBin1")
## fit the model
result <- hhh4(salmonella, model1)

## do sequential one-step-ahead predictions for the last 5 weeks
pred <- oneStepAhead(result, nrow(salmonella)-5, type="rolling",
                     which.start="final", verbose=FALSE)

\dontshow{
## test equivalence of parallelized version
if (.Platform$OS.type == "unix" && isTRUE(parallel::detectCores() > 1))
    stopifnot(identical(pred,
        oneStepAhead(result, nrow(salmonella)-5, type="rolling",
                     which.start="final", verbose=FALSE, cores=2)))
}

## oneStepAhead(..., type="final") just means fitted values
stopifnot(identical(
    unname(oneStepAhead(result, nrow(salmonella)-5,
                        type="final", verbose=FALSE)$pred),
    unname(tail(fitted(result), 5))))

## compute scores
(sc <- scores(pred, reverse = FALSE))

## scores with respect to the fitted values are similar
(scFitted <- scores(result, subset = nrow(salmonella)-c(4:0)))

\dontshow{
## test that scFitted is equivalent to scores(oneStepAhead(..., type = "final"))
stopifnot(all.equal(
    scFitted,
    scores(oneStepAhead(result, nrow(salmonella)-5, type="final", verbose=FALSE),
           reverse = FALSE),
    check.attributes = FALSE))
}

## plot a (non-randomized) PIT histogram for the predictions
with(pred, pit(x = observed, pdistr = "pnbinom", mu = pred, size = exp(psi)))


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

whichScores <- c("logs", "rps", "ses")
scores_i <- list()
meanScores <- NULL
for(i in seq(along.with=vals)){
  sc <- scores(get(vals[i]), which=whichScores, individual=TRUE)
  scores_i[[i]] <- sc
  meanScores <- rbind(meanScores,colMeans(sc, dims=2))
}

names(scores_i) <- nam
rownames(meanScores) <- nam

##comparison with best model B2 

compareWithBest <- function(best, whichModels, nPermut=9999, seed=1234){
  set.seed(seed)
  pVals <- NULL
  for(score in seq_along(whichScores)){
    p <- c()
    for(model in whichModels){
      if(model==best) p <- c(p,NA)
      else p <- c(p,permutationTest(scores_i[[model]][,,score],scores_i[[best]][,,score],
     plot=TRUE,nPermutation=nPermut, verbose=TRUE)$pVal.permut)
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
