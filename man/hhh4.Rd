\encoding{latin1}
\name{hhh4}
\alias{hhh4}

\title{Fitting HHH Models with Random Effects and Neighbourhood Structure}

\description{
  Fits an autoregressive Poisson or negative binomial model
  to a univariate or multivariate time series of counts.
  The characteristic feature of \code{hhh4} models is the additive
  decomposition of the conditional mean into \emph{epidemic} and
  \emph{endemic} components (Held et al, 2005).
  Log-linear predictors of covariates and random intercepts are allowed
  in all components; see the Details below.
  A general introduction to the \code{hhh4} modelling approach and its
  implementation is given in the \code{vignette("hhh4")}. Meyer et al
  (2017, Section 5, available as \code{vignette("hhh4_spacetime")})
  describe \code{hhh4} models for areal time series of infectious
  disease counts.
}

\usage{
hhh4(stsObj,
     control = list(
         ar = list(f = ~ -1, offset = 1, lag = 1),
         ne = list(f = ~ -1, offset = 1, lag = 1,
                   weights = neighbourhood(stsObj) == 1,
                   scale = NULL, normalize = FALSE),
         end = list(f = ~ 1, offset = 1),
         family = c("Poisson", "NegBin1", "NegBinM"),
         subset = 2:nrow(stsObj),
         optimizer = list(stop = list(tol=1e-5, niter=100),
                          regression = list(method="nlminb"),
                          variance = list(method="nlminb")),
         verbose = FALSE,
         start = list(fixed=NULL, random=NULL, sd.corr=NULL),
         data = list(t = stsObj@epoch - min(stsObj@epoch)),
         keep.terms = FALSE
     ),
     check.analyticals = FALSE)
}

\arguments{
  \item{stsObj}{object of class \code{"\linkS4class{sts}"} containing the (multivariate)
    count data time series.}
  \item{control}{a list containing the model specification and control arguments:
    \describe{
      \item{\code{ar}}{Model for the autoregressive component given as
        list with the following components: 
        \describe{
	  \item{f = ~ -1}{a formula specifying \eqn{\log(\lambda_{it})}{log(\lambda_it)}}
	  \item{offset = 1}{optional multiplicative offset, either 1 or
	    a matrix of the same dimension as \code{observed(stsObj)}}
	  \item{lag = 1}{a positive integer meaning autoregression on
	    \eqn{y_{i,t-lag}}}
	}
      }
      \item{\code{ne}}{Model for the neighbour-driven component given as
       list with the following components:
	\describe{
	  \item{f = ~ -1}{a formula specifying \eqn{\log(\phi_{it})}{log(\phi_it)}}
	  \item{offset = 1}{optional multiplicative offset, either 1 or
	    a matrix of the same dimension as \code{observed(stsObj)}}
	  \item{lag = 1}{a non-negative integer meaning dependency on
	    \eqn{y_{j,t-lag}}}
	  \item{weights = neighbourhood(stsObj) == 1}{
	    neighbourhood weights \eqn{w_{ji}}{w_ji}. The default
	    corresponds to the original formulation by Held et al
	    (2005), i.e., the spatio-temporal component incorporates an
	    unweighted sum over the lagged cases of the first-order
	    neighbours. See Paul et al (2008) and Meyer and Held (2014)
	    for alternative specifications, e.g.,
	    \code{\link{W_powerlaw}}.
	    Time-varying weights are possible by specifying an
	    array of \code{dim()} \code{c(nUnits, nUnits, nTime)}, where
	    \code{nUnits=ncol(stsObj)} and \code{nTime=nrow(stsObj)}.}
	  \item{scale = NULL}{
	    optional matrix of the same dimensions as \code{weights} (or
	    a vector of length \code{ncol(stsObj)}) to scale the
	    \code{weights} to \code{scale * weights}.
	  }
	  \item{normalize = FALSE}{
	    logical indicating if the (scaled) \code{weights} should be
	    normalized such that each row sums to 1.
	  }
	}
      }
      \item{\code{end}}{Model for the endemic component given as list
	with the following components
	\describe{
	  \item{f = ~ 1}{a formula specifying \eqn{\log(\nu_{it})}{log(\nu_it)}}
	  \item{offset = 1}{optional multiplicative offset \eqn{e_{it}}{e_it},
	    either 1 or a matrix of the same dimension as \code{observed(stsObj)}}
	}
      }
      \item{\code{family}}{Distributional family -- either \code{"Poisson"},
	or the Negative Binomial distribution. For the latter, the
	overdispersion parameter can be assumed to be the same for all
        units (\code{"NegBin1"}), to vary freely over all units
	(\code{"NegBinM"}), or to be shared by some units (specified by
        a factor of length \code{ncol(stsObj)} such that its number of
        levels determines the number of overdispersion parameters).
        Note that \code{"NegBinM"} is equivalent to
        \code{factor(colnames(stsObj), levels = colnames(stsObj))}.
      }
      \item{\code{subset}}{Typically \code{2:nrow(obs)} if model contains
	autoregression}
      \item{\code{optimizer}}{a list of three lists of control arguments.
	
	The \code{"stop"} list specifies two criteria for the outer
	optimization of regression and variance parameters: the relative
	\code{tol}erance for parameter change using the criterion 
	\code{max(abs(x[i+1]-x[i])) / max(abs(x[i]))},
	and the maximum number \code{niter} of outer iterations.
	
	Control arguments for the single optimizers are specified in the
	lists named \code{"regression"} and \code{"variance"}.
	\code{method="nlminb"} is the default optimizer for both (taking
	advantage of the analytical Fisher information matrices), however,
	the \code{method}s from \code{\link{optim}} may also be specified
	(as well as \code{"\link{nlm}"} but that one is not recommended here).
	Especially for the variance updates, Nelder-Mead optimization
	(\code{method="Nelder-Mead"}) is an attractive alternative.
	All other elements of these two lists are passed as
	\code{control} arguments to the chosen \code{method}, e.g., if
	\code{method="nlminb"} adding \code{iter.max=50} increases the 
	maximum number of inner iterations from 20 (default) to 50.
      }
      \item{\code{verbose}}{non-negative integer (usually in the range
	\code{0:3}) specifying the amount of tracing information to be
	output during optimization.}
      \item{\code{start}}{a list of initial parameter values replacing
	initial values set via \code{\link{fe}} and \code{\link{ri}}.
	Since \pkg{surveillance} 1.8-2, named vectors are matched
	against the coefficient names in the model (where unmatched
	start values are silently ignored), and need not be complete,
	e.g., \code{start = list(fixed = c("-log(overdisp)" = 0.5))}
	(default: 2) for a \code{family = "NegBin1"} model.
	In contrast, an unnamed start vector must specify the full set
	of parameters as used by the model.}
      \item{\code{data}}{a named list of covariates that are to be
	included as fixed effects (see \code{\link{fe}}) in any of the 3
	component formulae.
	By default, the time variable \code{t} is available and used for
	seasonal effects created by \code{\link{addSeason2formula}}.
	In general, covariates in this list can be either vectors of
	length \code{nrow(stsObj)} interpreted as time-varying but
	common across all units, or matrices of the same dimension as
	the disease counts \code{observed(stsObj)}.}
      \item{\code{keep.terms}}{logical indicating if the terms object
	used in the fit is to be kept as part of the returned object.
	This is usually not necessary, since the terms object is
	reconstructed by the \code{\link{terms}}-method for class
	\code{"hhh4"} if necessary (based on \code{stsObj} and
	\code{control}, which are both part of the returned
	\code{"hhh4"} object).}
    }
    The auxiliary function \code{\link{makeControl}} might be useful to
    create such a list of control parameters.
  }
  \item{check.analyticals}{logical (or a subset of
    \code{c("numDeriv", "maxLik")}), indicating if (how) the implemented
    analytical score vector and Fisher information matrix should be
    checked against numerical derivatives at the parameter starting values,
    using the packages \pkg{numDeriv} and/or \pkg{maxLik}. If activated,
    \code{hhh4} will return a list containing the analytical and numerical
    derivatives for comparison (no ML estimation will be performed). 
    This is mainly intended for internal use by the package developers.}
}

\value{
  \code{hhh4} returns an object of class \code{"hhh4"},
  which is a list containing the following components:
  \item{coefficients}{named vector with estimated (regression) parameters of the model}
  \item{se}{estimated standard errors (for regression parameters)}
  \item{cov}{covariance matrix (for regression parameters)}
  \item{Sigma}{estimated variance-covariance matrix of random effects}
  \item{Sigma.orig}{estimated variance parameters on internal scale used
    for optimization}
  \item{Sigma.cov}{inverse of marginal Fisher information (on internal
    scale), i.e., the asymptotic covariance matrix of \code{Sigma.orig}}
  \item{call}{ the matched call }
  \item{dim}{ vector with number of fixed and random effects in the model }
  \item{loglikelihood}{(penalized) loglikelihood evaluated at the MLE}
  \item{margll}{ (approximate) log marginal likelihood should the model contain random effects  }
  \item{convergence}{logical. Did optimizer converge?}
  \item{fitted.values}{fitted mean values \eqn{\mu_{i,t}}{\mu_it}}
  \item{control}{control object of the fit}
  \item{terms}{the terms object used in the fit if \code{keep.terms = TRUE}
    and \code{NULL} otherwise}
  \item{stsObj}{ the supplied \code{stsObj} }
  \item{lags}{named integer vector of length two containing the lags
    used for the epidemic components \code{"ar"} and \code{"ne"},
    respectively. The corresponding lag is \code{NA} if the component
    was not included in the model.}
  \item{nObs}{number of observations used for fitting the model}
  \item{nTime}{ number of time points used for fitting the model }
  \item{nUnit}{ number of units (e.g. areas) used for fitting the model}
  \item{runtime}{the \code{\link{proc.time}}-queried time taken
    to fit the model, i.e., a named numeric vector of length 5 of class
    \code{"proc_time"}}
}

\details{
  An endemic-epidemic multivariate time-series model for infectious
  disease counts \eqn{Y_{it}}{Y_it} from units \eqn{i=1,\dots,I} during
  periods \eqn{t=1,\dots,T} was proposed by Held et al (2005) and was
  later extended in a series of papers (Paul et al, 2008; Paul and Held,
  2011; Held and Paul, 2012; Meyer and Held, 2014).
  In its most general formulation, this so-called \code{hhh4} (or HHH or
  \eqn{H^3} or triple-H) model assumes that, conditional on past
  observations, \eqn{Y_{it}}{Y_it} has a Poisson or negative binomial
  distribution with mean
  \deqn{\mu_{it} = \lambda_{it} y_{i,t-1} + 
                   \phi_{it} \sum_{j\neq i} w_{ji} y_{j,t-1} +
                   e_{it} \nu_{it}  }{%
        \mu_it = \lambda_it y_i,t-1 + 
                 \phi_it sum_(j != i) w_ji y_j,t-1 +
                   e_it \nu_it  }
  In the case of a negative binomial model, the conditional 
  variance is \eqn{\mu_{it}(1+\psi_i\mu_{it})}{\mu_it(1+\psi_i*\mu_it)} 
  with overdispersion parameters \eqn{\psi_i > 0} (possibly shared
  across different units, e.g., \eqn{\psi_i\equiv\psi}{\psi_i=\psi}).
  Univariate time series of counts \eqn{Y_t} are supported as well, in
  which case \code{hhh4} can be regarded as an extension of
  \code{\link[MASS]{glm.nb}} to account for autoregression.
  See the Examples below for a comparison of an endemic-only
  \code{hhh4} model with a corresponding \code{glm.nb}.
  
  The three unknown quantities of the mean \eqn{\mu_{it}}{\mu_it},
  \itemize{
  \item \eqn{\lambda_{it}}{\lambda_it} in the autoregressive (\code{ar}) component, 
  \item \eqn{\phi_{it}}{\phi_it} in the neighbour-driven (\code{ne}) component, and
  \item \eqn{\nu_{it}}{\nu_it} in the endemic (\code{end}) component,
  }
  are log-linear predictors incorporating time-/unit-specific
  covariates. They may also contain unit-specific random intercepts
  as proposed by Paul and Held (2011). The endemic mean is usually
  modelled proportional to a unit-specific offset \eqn{e_{it}}{e_it}
  (e.g., population numbers or fractions); it is possible to include
  such multiplicative offsets in the epidemic components as well.
  The \eqn{w_{ji}}{w_ji} are transmission weights reflecting the flow of
  infections from unit \eqn{j} to unit \eqn{i}.
  In spatial \code{hhh4} applications, the \dQuote{units} refer to
  geographical regions and the weights could be derived from movement
  network data. Alternatively, the weights can be estimated
  parametrically as a function of adjacency order (Meyer and Held, 2014).

  (Penalized) Likelihood inference for such \code{hhh4} models has been
  established by Paul and Held (2011) with extensions for parametric
  neighbourhood weights by Meyer and Held (2014).
  Supplied with the analytical score function and Fisher information,
  the function \code{hhh4} by default uses the quasi-Newton algorithm
  available through \code{\link{nlminb}} to maximize the log-likelihood.
  Convergence is usually fast even for a large number of parameters.
  If the model contains random effects, the penalized and marginal
  log-likelihoods are maximized alternately until convergence.
}

\seealso{
  See the special functions \code{\link{fe}}, \code{\link{ri}} and the
  examples below for how to specify unit-specific effects.

  Further details on the modelling approach and illustrations of its
  implementation can be found in \code{vignette("hhh4")} and
  \code{vignette("hhh4_spacetime")}.
}

\author{Michaela Paul, Sebastian Meyer, Leonhard Held}

\examples{
######################
## Univariate examples
######################

### weekly counts of salmonella agona cases, UK, 1990-1995

data("salmonella.agona")
## convert old "disProg" to new "sts" data class
salmonella <- disProg2sts(salmonella.agona)
salmonella
plot(salmonella)

## generate formula for an (endemic) time trend and seasonality
f.end <- addSeason2formula(f = ~1 + t, S = 1, period = 52)
f.end
## specify a simple autoregressive negative binomial model
model1 <- list(ar = list(f = ~1), end = list(f = f.end), family = "NegBin1")
## fit this model to the data
res <- hhh4(salmonella, model1)
## summarize the model fit
summary(res, idx2Exp=1, amplitudeShift=TRUE, maxEV=TRUE)
plot(res)
plot(res, type = "season", components = "end")


### weekly counts of meningococcal infections, Germany, 2001-2006

data("influMen")
fluMen <- disProg2sts(influMen)
meningo <- fluMen[, "meningococcus"]
meningo
plot(meningo)

## again a simple autoregressive NegBin model with endemic seasonality
meningoFit <- hhh4(stsObj = meningo, control = list(
    ar = list(f = ~1),
    end = list(f = addSeason2formula(f = ~1, S = 1, period = 52)),
    family = "NegBin1"
))

summary(meningoFit, idx2Exp=TRUE, amplitudeShift=TRUE, maxEV=TRUE)
plot(meningoFit)
plot(meningoFit, type = "season", components = "end")


########################
## Multivariate examples
########################

### bivariate analysis of influenza and meningococcal infections
### (see Paul et al, 2008)

plot(fluMen, same.scale = FALSE)
     
## Fit a negative binomial model with
## - autoregressive component: disease-specific intercepts
## - neighbour-driven component: only transmission from flu to men
## - endemic component: S=3 and S=1 sine/cosine pairs for flu and men, respectively
## - disease-specific overdispersion

WfluMen <- neighbourhood(fluMen)
WfluMen["meningococcus","influenza"] <- 0
WfluMen
f.end_fluMen <- addSeason2formula(f = ~ -1 + fe(1, which = c(TRUE, TRUE)),
                                  S = c(3, 1), period = 52)
f.end_fluMen
fluMenFit <- hhh4(fluMen, control = list(
    ar = list(f = ~ -1 + fe(1, unitSpecific = TRUE)),
    ne = list(f = ~ 1, weights = WfluMen),
    end = list(f = f.end_fluMen),
    family = "NegBinM"))
summary(fluMenFit, idx2Exp=1:3)
plot(fluMenFit, type = "season", components = "end", unit = 1)
plot(fluMenFit, type = "season", components = "end", unit = 2)
\dontshow{
    ## regression test for amplitude/shift transformation of sine-cosine pairs
    ## coefficients were wrongly matched in surveillance < 1.18.0
    stopifnot(coef(fluMenFit, amplitudeShift = TRUE)["end.A(2 * pi * t/52).meningococcus"] == sqrt(sum(coef(fluMenFit)[paste0("end.", c("sin","cos"), "(2 * pi * t/52).meningococcus")]^2)))
}


### weekly counts of measles, Weser-Ems region of Lower Saxony, Germany

data("measlesWeserEms")
measlesWeserEms
plot(measlesWeserEms)  # note the two districts with zero cases

## we could fit the same simple model as for the salmonella cases above
model1 <- list(
    ar = list(f = ~1),
    end = list(f = addSeason2formula(~1 + t, period = 52)),
    family = "NegBin1"
)
measlesFit <- hhh4(measlesWeserEms, model1)
summary(measlesFit, idx2Exp=TRUE, amplitudeShift=TRUE, maxEV=TRUE)

## but we should probably at least use a population offset in the endemic
## component to reflect heterogeneous incidence levels of the districts,
## and account for spatial dependence (here just using first-order adjacency)
measlesFit2 <- update(measlesFit,
    end = list(offset = population(measlesWeserEms)),
    ne = list(f = ~1, weights = neighbourhood(measlesWeserEms) == 1))
summary(measlesFit2, idx2Exp=TRUE, amplitudeShift=TRUE, maxEV=TRUE)
plot(measlesFit2, units = NULL, hide0s = TRUE)

## 'measlesFit2' corresponds to the 'measlesFit_basic' model in
## vignette("hhh4_spacetime"). See there for further analyses,
## including vaccination coverage as a covariate,
## spatial power-law weights, and random intercepts.


\dontrun{
### last but not least, a more sophisticated (and time-consuming)
### analysis of weekly counts of influenza from 140 districts in
### Southern Germany (originally analysed by Paul and Held, 2011,
### and revisited by Held and Paul, 2012, and Meyer and Held, 2014)

data("fluBYBW")
plot(fluBYBW, type = observed ~ time)
plot(fluBYBW, type = observed ~ unit,
     ## mean yearly incidence per 100.000 inhabitants (8 years)
     population = fluBYBW@map$X31_12_01 / 100000 * 8)

## For the full set of models for data("fluBYBW") as analysed by
## Paul and Held (2011), including predictive model assessement
## using proper scoring rules, see the (computer-intensive)
## demo("fluBYBW") script:
demoscript <- system.file(file.path("demo", "fluBYBW.R"),
                          package = "surveillance")
demoscript
#file.show(demoscript)

## Here we fit the improved power-law model of Meyer and Held (2014)
## - autoregressive component: random intercepts + S = 1 sine/cosine pair
## - neighbour-driven component: random intercepts + S = 1 sine/cosine pair
##   + population gravity with normalized power-law weights
## - endemic component: random intercepts + trend + S = 3 sine/cosine pairs
## - random intercepts are iid but correlated between components
f.S1 <- addSeason2formula(
    ~-1 + ri(type="iid", corr="all"),
    S = 1, period = 52)
f.end.S3 <- addSeason2formula(
    ~-1 + ri(type="iid", corr="all") + I((t-208)/100),
    S = 3, period = 52)

## for power-law weights, we need adjaceny orders, which can be
## computed from the binary adjacency indicator matrix
nbOrder1 <- neighbourhood(fluBYBW)
neighbourhood(fluBYBW) <- nbOrder(nbOrder1, 15)

## full model specification
fluModel <- list(
    ar = list(f = f.S1),
    ne = list(f = update.formula(f.S1, ~ . + log(pop)),
              weights = W_powerlaw(maxlag=max(neighbourhood(fluBYBW)),
                                   normalize = TRUE, log = TRUE)),
    end = list(f = f.end.S3, offset = population(fluBYBW)),
    family = "NegBin1", data = list(pop = population(fluBYBW)),
    optimizer = list(variance = list(method = "Nelder-Mead")),
    verbose = TRUE)

## CAVE: random effects considerably increase the runtime of model estimation
## (It is usually advantageous to first fit a model with simple intercepts
## to obtain reasonable start values for the other parameters.)
set.seed(1)  # because random intercepts are initialized randomly
fluFit <- hhh4(fluBYBW, fluModel)

summary(fluFit, idx2Exp = TRUE, amplitudeShift = TRUE)

plot(fluFit, type = "fitted", total = TRUE)

plot(fluFit, type = "season")
range(plot(fluFit, type = "maxEV"))

plot(fluFit, type = "maps", prop = TRUE)

gridExtra::grid.arrange(
    grobs = lapply(c("ar", "ne", "end"), function (comp)
        plot(fluFit, type = "ri", component = comp, main = comp,
             exp = TRUE, sub = "multiplicative effect")),
    nrow = 1, ncol = 3)

plot(fluFit, type = "neweights", xlab = "adjacency order")
}


########################################################################
## An endemic-only "hhh4" model can also be estimated using MASS::glm.nb
########################################################################

## weekly counts of measles, Weser-Ems region of Lower Saxony, Germany
data("measlesWeserEms")

## fit an endemic-only "hhh4" model
## with time covariates and a district-specific offset
hhh4fit <- hhh4(measlesWeserEms, control = list(
    end = list(f = addSeason2formula(~1 + t, period = measlesWeserEms@freq),
               offset = population(measlesWeserEms)),
    ar = list(f = ~-1), ne = list(f = ~-1), family = "NegBin1",
    subset = 1:nrow(measlesWeserEms)
))
summary(hhh4fit)

## fit the same model using MASS::glm.nb
measlesWeserEmsData <- as.data.frame(measlesWeserEms, tidy = TRUE)
measlesWeserEmsData$t <- c(hhh4fit$control$data$t)
glmnbfit <- MASS::glm.nb(
    update(formula(hhh4fit)$end, observed ~ . + offset(log(population))),
    data = measlesWeserEmsData
)
summary(glmnbfit)

## Note that the overdispersion parameter is parametrized inversely.
## The likelihood and point estimates are all the same.
## However, the variance estimates are different: in glm.nb, the parameters
## are estimated conditional on the overdispersion theta.

\dontshow{
stopifnot(
    all.equal(logLik(hhh4fit), logLik(glmnbfit)),
    all.equal(1/coef(hhh4fit)[["overdisp"]], glmnbfit$theta, tolerance = 1e-6),
    all.equal(coef(hhh4fit)[1:4], coef(glmnbfit),
              tolerance = 1e-6, check.attributes = FALSE),
    all.equal(c(residuals(hhh4fit)), residuals(glmnbfit),
              tolerance = 1e-6, check.attributes = FALSE)
)
}
}

\references{
  Held, L., \enc{Höhle}{Hoehle}, M. and Hofmann, M. (2005):
    A statistical framework for the analysis of multivariate infectious
    disease surveillance counts.
    \emph{Statistical Modelling}, \bold{5} (3), 187-199.
    \doi{10.1191/1471082X05st098oa}
    
  Paul, M., Held, L. and Toschke, A. M. (2008):
    Multivariate modelling of infectious disease surveillance data.
    \emph{Statistics in Medicine}, \bold{27} (29), 6250-6267.
    \doi{10.1002/sim.4177}

  Paul, M. and Held, L. (2011):
    Predictive assessment of a non-linear random effects model for
    multivariate time series of infectious disease counts.
    \emph{Statistics in Medicine}, \bold{30} (10), 1118-1136.
    \doi{10.1002/sim.4177}

  Held, L. and Paul, M. (2012):
    Modeling seasonality in space-time infectious disease surveillance data.
    \emph{Biometrical Journal}, \bold{54} (6), 824-843.
    \doi{10.1002/bimj.201200037}
    
  Meyer, S. and Held, L. (2014):
    Power-law models for infectious disease spread.
    \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.
    \doi{10.1214/14-AOAS743}

  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
    Spatio-temporal analysis of epidemic phenomena using the \R package
    \pkg{surveillance}.
    \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
    \doi{10.18637/jss.v077.i11}
}

\keyword{ts}
\keyword{regression}
