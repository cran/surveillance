\encoding{latin1}
\name{twinSIR_simulation}
\alias{simEpidata}
\alias{simulate.twinSIR}

\title{
Simulation of Epidemic Data
}

\description{
  This function simulates the infection (and removal) times of 
  an epidemic.  Besides the classical SIR type of epidemic, also SI, SIRS and
  SIS epidemics are supported.  Simulation works via the conditional intensity 
  of infection of an individual, given some (time varying) endemic covariates 
  and/or some distance functions (epidemic components) as well as the fixed 
  positions of the individuals.  The lengths of the infectious and removed 
  periods are generated following a pre-specified function (can be
  deterministic).
  
  The \code{\link{simulate}} method for objects of class
  \code{"\link{twinSIR}"} simulates new epidemic data using the model and
  the parameter estimates of the fitted object.
}

\usage{
simEpidata(formula, data, id.col, I0.col, coords.cols, subset,
           beta, h0, f = list(), w = list(), alpha, infPeriod,
           remPeriod = function(ids) rep(Inf, length(ids)),
           end = Inf, trace = FALSE, .allocate = NULL)

\method{simulate}{twinSIR}(object, nsim = 1, seed = 1,
         infPeriod = NULL, remPeriod = NULL,
         end = diff(range(object$intervals)), trace = FALSE, .allocate = NULL,
         data = object$data, ...)
}

\arguments{
  \item{formula}{
    an object of class \code{"\link{formula}"} (or one that can be coerced to
    that class): a symbolic description of the intensity model to be estimated.  
    The details of model specification are given under Details.
  }
  \item{data}{
    a data.frame containing the variables in \code{formula} and the variables
    specified by \code{id.col}, \code{I0.col} and \code{coords.col} (see below).
    It represents the \dQuote{history} of the endemic covariates to use for the
    simulation.  The form is similar to and can be an object of class
    \code{"\link{epidata}"}.  The simulation period is split up into
    \emph{consecutive} intervals of constant endemic covariables.  The 
    data frame consists of a block of N (number of individuals) rows for each of
    those time intervals (all rows in a block share the same start and stop 
    values... therefore the name \dQuote{block}), where there is one row per 
    individual in the block.  Each row describes the (fixed) state of the 
    endemic covariates of the individual during the time interval given by the 
    start and stop columns (specified through the lhs of \code{formula}).
    
    For the \code{simulate} method of class \code{"twinSIR"} this should be 
    the object of class \code{"\link{epidata}"} used for the fit.  This is a
    part of the return value of the function \code{twinSIR}, if called with 
    argument \code{keep.data} set to \code{TRUE}.
  }
  \item{id.col}{
    only if \code{data} does not inherit from \code{epidata}:
    single index of the \code{id} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).\cr
    The \code{id} column identifies the individuals in the data-frame.  It will
    be converted to a factor variable and its levels serve also to identify 
    individuals as argument to the \code{infPeriod} function.
  }
  \item{I0.col}{
    only if \code{data} does not inherit from \code{epidata}:
    single index of the \code{I0} column in \code{data}.  Can be numeric
    (by column number), character (by column name) or \code{NULL}.\cr
    The \code{I0} column indicates if an individual is initially infectious,
    i.e. it is already infectious at the beginning of the first time block.
    Setting \code{I0.col = NULL} is short for \dQuote{there are no initially
    infectious individuals}. Otherwise, the variable must be logical or in
    0/1-coding.  As this variable is constant over time the initially
    infectious individuals are derived from the first time block only.
}
  \item{coords.cols}{
    only if \code{data} does not inherit from \code{epidata}:
    index\emph{es} of the \code{coords} column\emph{s} in \code{data}.  Can be a
    numeric (by column number), a character (by column name) vector or
    \code{NULL}.\cr
    These columns contain the coordinates of the individuals.  It must be
    emphasized that the functions in this package currently assume \emph{fixed
    positions} of the individuals during the whole epidemic.  Thus, an
    individual has the same coordinates in every block.  For simplicity, the
    coordinates are derived from the first time block only.  The epidemic
    covariates are calculated based on the Euclidian distance between the
    individuals, see \code{f}.
  }
  \item{subset}{
    an optional vector specifying a subset of the covariate history to be used 
    in the simulation.
  }
  \item{beta}{
    numeric vector of length equal the number of endemic (\code{cox}) terms on
    the rhs of \code{formula}.  It contains the effects of the endemic predictor
    (excluding the log-baseline \code{h0}, see below) in the same order as in 
    the formula.
  }
  \item{h0}{
    \emph{either} a single number to specify a constant baseline hazard
    (equal to \code{exp(h0)}) \emph{or} a list of functions named
    \code{exact} and \code{upper}.  In the latter case, \code{h0$exact}
    is the true log-baseline hazard function and \code{h0$upper} is a
    \emph{piecewise constant upper bound} for \code{h0$exact}.  The function
    \code{h0$upper} must inherit from \code{\link{stepfun}} with
    \code{right=FALSE}.  Theoretically, the intensity function is
    left-continuous, thus \code{right=TRUE} would be adequate, but in the
    implementation, when we evaluate the intensity at the
    \code{\link{knots}} (change points) of \code{h0$upper} we need its value
    for the subsequent interval.
  }
  \item{f, w}{
    see \code{\link{as.epidata}}.
  }
  \item{alpha}{
    a named numeric vector of coefficients for the epidemic
    covariates generated by \code{f} and \code{w}. The names are matched
    against \code{names(f)} and \code{names(w)}.
    Remember that \code{alpha >= 0}.
  }
  \item{infPeriod}{
    a function generating lengths of infectious periods.  It should take one
    parameter (e.g. \code{ids}), which is a character vector of id's of 
    individuals, and return appropriate infection periods for those
    individuals.  Therefore, the value of the function should be of length
    \code{length(ids)}.  For example, for independent and identically
    distributed infection periods following \eqn{Exp(1)},
    the generating function is \code{function(ids) rexp(length(ids), rate=1)}.  
    For a constant infectious period of length c, it is sufficient to set
    \code{function (x) {c}}.\cr
    For the \code{simulate} method of class \code{"twinSIR"} only, this can
    also be \code{NULL} (the default), which means that the observed infectious
    periods of infected individuals are re-used when simulating a new epidemic
    and individuals with missing infectious periods (i.e. infection and 
    recovery was not observed) are attributed to the mean observed
    infectious period.
    
    Note that it is even possible to simulate an SI-epidemic by setting

    \code{infPeriod = function (x) {Inf}}

    In other words: once an individual became
    infected it spreads the disease forever, i.e. it will never be removed.
  }
  \item{remPeriod}{
    a function generating lengths of removal periods.  Per default, once an
    individual was removed it will stay in this state forever (\code{Inf}).
    Therefore, it will not become at-risk (S) again and re-infections are not
    possible.  Alternatively, always returning 0 as length of the removal 
    period corresponds to a SIS epidemic.  Any other values correspond to SIRS.
    Note that \code{end} should be set to a finite value in these cases.
  }
  \item{end}{
    a single positive numeric value specifying the time point at which the
    simulation should be forced to end.  By default, this is \code{Inf}, i.e.
    the simulation continues until there is no susceptible individual left.\cr
    For the \code{simulate} method of class \code{"twinSIR"} the default is
    to have equal simulation and observation periods.
  }
  \item{trace}{
    logical (or integer) indicating if (or how often) the sets of susceptible 
    and infected individuals as well as the rejection indicator (of the 
    rejection sampling step) should be \code{cat}ed.  Defaults to \code{FALSE}.
  }
  \item{.allocate}{
    number of blocks to initially allocate for the event history (i.e.
    \code{.allocate*N} rows). By default (\code{NULL}), this number is set to
    \code{max(500, ceiling(nBlocks/100)*100)}, i.e. 500 but at least the
    number of blocks in \code{data} (rounded to the next multiple of 100).
    Each time the simulated epidemic exceeds the allocated space, the
    event history will be enlarged by \code{.allocate} blocks.
  }
  \item{object}{
    an object of class \code{"twinSIR"}.  This must contain the original
    \code{data} used for the fit (see \code{data}).
  }
  \item{nsim}{
    number of epidemics to simulate.  Defaults to 1.
  }
  \item{seed}{
    an integer that will be used in the call to \code{\link{set.seed}} before
    simulating the epidemics.
  }
  \item{\dots}{
    unused (argument of the generic).
  }
}

\details{
  A model is specified through the \code{formula}, which has the form

  \code{cbind(start, stop) ~ cox(endemicVar1) * cox(endemicVar2)},
  
  i.e. the right hand side has the usual form as in \code{\link{lm}}, but
  all variables are marked as being endemic by the special function
  \code{\link{cox}}.  The effects of those predictor terms are specified by
  \code{beta}.  The left hand side of the formula denotes the start
  and stop columns in \code{data}.  This can be omitted, if \code{data} inherits
  from class \code{"epidata"} in which case \code{cbind(start, stop)} will be
  used.  The epidemic model component is specified by the arguments
  \code{f} and \code{w} (and the associated coefficients \code{alpha}).

  If the epidemic model component is empty and \code{infPeriod}
  always returns \code{Inf}, then one actually simulates from a pure Cox model.

  The simulation algorithm used is \emph{Ogata's modified thinning}.
  For details, see \enc{Höhle}{Hoehle} (2009), Section 4.
}

\value{
  An object of class \code{"simEpidata"}, which is a \code{data.frame} with the
  columns \code{"id"}, \code{"start"}, \code{"stop"}, \code{"atRiskY"},
  \code{"event"}, \code{"Revent"} and the coordinate columns (with the original
  names from \code{data}), which are all obligatory.  These columns are followed
  by all the variables appearing on the rhs of the \code{formula}.  Last but not
  least, the generated columns with epidemic covariates corresponding to the
  functions in the lists \code{f} and \code{w} are appended.

  Note that objects of class \code{"simEpidata"} also inherit from class
  \code{"\link{epidata}"}, thus all \code{"\link{epidata}"} methods can be
  applied.
  
  The \code{data.frame} is given the additional \emph{attributes}
  \item{"eventTimes"}{
    numeric vector of infection time points (sorted chronologically).
  }
  \item{"timeRange"}{
    numeric vector of length 2: \code{c(min(start), max(stop))}.
  }
  \item{"coords.cols"}{
    numeric vector containing the column indices of the coordinate columns in
    the resulting data-frame.
  }
  \item{"f"}{
    this equals the argument \code{f}.
  }
  \item{"w"}{
    this equals the argument \code{w}.
  }
  \item{"config"}{
    a list with elements \code{h0 = h0$exact}, \code{beta} and \code{alpha}.
  }
  \item{call}{the matched call.}
  \item{terms}{the \code{terms} object used.}

  If \code{nsim > 1} epidemics are simulated by the
  \code{simulate}-method for fitted \code{"twinSIR"} models, these are
  returned in a list.
}

\references{
  \enc{Höhle}{Hoehle}, M. (2009),
  Additive-Multiplicative Regression Models for Spatio-Temporal
  Epidemics, Biometrical Journal, 51(6):961-978.
}

\author{ 
Sebastian Meyer and Michael \enc{Höhle}{Hoehle}
}

\seealso{
The \code{\link{plot.epidata}} and \code{\link{animate.epidata}} methods
for plotting and animating (simulated) epidemic data, respectively.
The \code{\link{intensityplot.simEpidata}} method for plotting paths of
infection intensities.

Function \code{\link{twinSIR}} for fitting spatio-temporal epidemic intensity
models to epidemic data.
}

\examples{
## Generate a data frame containing a hypothetic population with 100 individuals
set.seed(1234)
n <- 100
pos <- matrix(rnorm(n*2), ncol=2, dimnames=list(NULL, c("x", "y")))
pop <- data.frame(id=1:n, x=pos[,1], y=pos[,2], 
                  gender=sample(0:1, n, replace=TRUE),
                  I0col=c(rep(1,3),rep(0,n-3)), # 3 initially infectious
                  start=rep(0,n), stop=rep(Inf,n))

## Simulate an SIR epidemic in this population
set.seed(123)
infPeriods <- setNames(c(1:3/10, rexp(n-3, rate=1)), 1:n)
epi <- simEpidata(
    cbind(start,stop) ~ cox(gender), data = pop,
    id = "id", I0.col = "I0col", coords.cols = c("x","y"),
    beta = c(-2), h0 = -1, alpha = c(B1=0.1), f = list(B1=function(u) u<=1),
    infPeriod = function(ids) infPeriods[ids],
    ##remPeriod = function(ids) rexp(length(ids), rate=0.1), end = 30   # -> SIRS
)

## extract event times by id
head(summary(epi)$byID)

## Plot the numbers of susceptible, infectious and removed individuals
plot(epi)


## load the 1861 Hagelloch measles epidemic
data("hagelloch")
summary(hagelloch)
plot(hagelloch)

## fit a simplistic twinSIR model
fit <- twinSIR(~ household, data = hagelloch)

## simulate a new epidemic from the above model
## with simulation period = observation period, re-using observed infPeriods
sim1 <- simulate(fit, data = hagelloch)
plot(sim1)

## check if we find similar parameters in the simulated epidemic
fitsim1 <- update(fit, data = sim1)
cbind(base = coef(fit), new = coef(fitsim1))

if (surveillance.options("allExamples"))
{
  ## simulate only 10 days, using random infPeriods ~ Exp(0.1)
  sim2 <- simulate(fit, data = hagelloch, seed = 2, end = 10,
    infPeriod = function(ids) rexp(length(ids), rate = 0.1))
  plot(sim2)

  ## simulate from a different model with manually specified parameters
  set.seed(321)
  simepi <- simEpidata(~ cox(AGE), data = hagelloch,
      beta = c(0.1), h0 = -4, alpha = c(household = 0.05),
      f = list(household = function(u) u == 0),
      infPeriod = function(ids) rexp(length(ids), rate=1/8))
  plot(simepi)
  intensityplot(simepi)

  ## see if we correctly estimate the parameters
  fitsimepi <- twinSIR(~ cox(AGE) + household, data = simepi)
  cbind(true = c(0.05, -4, 0.1), est = coef(fitsimepi), confint(fitsimepi))
}
}

\keyword{datagen}
\keyword{models}
