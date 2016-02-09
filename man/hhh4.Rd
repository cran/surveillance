\encoding{latin1}
\name{hhh4}
\alias{hhh4}

\title{Fitting HHH Models with Random Effects and Neighbourhood Structure}

\description{
  Fits a Poisson or negative binomial model with conditional mean
  \deqn{\mu_{it} = \lambda_{it} y_{i,t-1} + 
                   \phi_{it} \sum_{j\neq i} w_{ji} y_{j,t-1} +
                   e_{it} \nu_{it}  }{%
        \mu_it = \lambda_it y_i,t-1 + 
                 \phi_it sum_(j != i) w_ji y_j,t-1 +
                   e_it \nu_it  }
  containing epidemic and endemic components to a multivariate time series
  of counts \eqn{Y_{it}}{Y_it} (unit \eqn{i}, period \eqn{t}).
  Univariate count time series are supported as well.
  In the case of a negative binomial model, the conditional 
  variance is \eqn{\mu_{it}(1+\psi\mu_{it})}{\mu_it(1+\psi*\mu_it)} 
  with overdispersion parameter \eqn{\psi}. 
  The three unknown quantities of the mean \eqn{\mu_{it}}{\mu_it},
  \itemize{
  \item \eqn{\lambda_{it}}{\lambda_it} in the autoregressive (\code{ar}) component, 
  \item \eqn{\phi_{it}}{\phi_it} in the neighbor-driven (\code{ne}) component, and
  \item \eqn{\nu_{it}}{\nu_it} in the endemic (\code{end}) component,
  }
  are log-linear predictors incorporating time-/unit-specific
  covariates. They may contain unit-specific random intercepts
  (\code{\link{ri}}) as proposed by Paul and Held (2011).
  The \eqn{e_{it}}{e_it} is a (multiplicative) endemic offset; it is
  also possible to include such offsets in the epidemic components.
  The \eqn{w_{ji}}{w_ji} are neighbourhood weights, which can either be
  prespecified or estimated parametrically as proposed by Meyer and Held (2014).
  The implementation is illustrated in the \code{vignette("hhh4")} and
  in Meyer et al. (2016, Section 5).
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
      \item{\code{ne}}{Model for the neighbor-driven component given as
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
  For further details see \code{vignette("hhh4")} and the references.
}

\seealso{
  \code{\link{algo.hhh}}, 
  \code{\link{fe}}, \code{\link{ri}}
}

\author{M. Paul, S. Meyer, and L. Held}

\examples{
#####################################################################
# Fit some models from ?algo.hhh
#####################################################################

## univariate salmonella agona data
data(salmonella.agona)
# convert to sts class
salmonella <- disProg2sts(salmonella.agona)

# generate formula for temporal and seasonal trends
f.end <- addSeason2formula(f = ~ 1 + t, S=1, period=52)
model1 <- list(ar = list(f = ~ 1), end = list(f =f.end),
               family = "NegBin1")
# run model
res <- hhh4(salmonella, model1)
summary(res, idx2Exp=1, amplitudeShift=TRUE)


## multivariate time series: 
# measles cases in Lower Saxony, Germany
data(measles.weser)
measles <- disProg2sts(measles.weser)

# same model as above
summary(hhh4(measles, control=model1))

# now use region-specific intercepts in endemic component
f.end2 <- addSeason2formula(f = ~ -1 + fe(1, unitSpecific = TRUE) + t,
                            S = 1, period = 52)
model2 <- list(ar = list(f = ~ 1), 
               end = list(f = f.end2, offset = population(measles)),
               family = "NegBin1")
# run model
summary(hhh4(measles, control=model2), idx2Exp=1, amplitudeShift=TRUE)

# include autoregressive parameter phi for adjacent "Kreise"
# no linear trend in endemic component
f.end3 <- addSeason2formula(f = ~ -1 + fe(1, unitSpecific = TRUE), 
                            S = 1, period = 52)
model3 <- list(ar = list(f = ~ 1),
               ne = list(f = ~ 1),
               end = list(f = f.end3, offset = population(measles)),
               family = "NegBin1")
# run model
res3 <- hhh4(measles, control=model3)
summary(res3, idx2Exp=1:2, amplitudeShift=TRUE)

\dontshow{
## check that neighbourhood weights array yields same result
.neweights <- array(neighbourhood(measles),
                   dim = c(rep(ncol(measles),2),nrow(measles)),
                   dimnames = c(dimnames(neighbourhood(measles)), list(NULL)))
res3_tv <- hhh4(measles, control = modifyList(model3,
                list(ne=list(weights=.neweights))))
stopifnot(all.equal(res3, res3_tv, ignore = "control"))
}

\dontrun{
######################################################################
# Fit the models from the Paul & Held (2011) paper for the influenza data
# from Bavaria and Baden-Wuerttemberg (this takes some time!)
# For further documentation see also the vignette.
######################################################################

data("fluBYBW") 

###############################################################
## generate formula for temporal and seasonal trends
f.end <- addSeason2formula(f = ~ -1 + ri(type="iid", corr="all") + 
                               I((t-208)/100), S=3, period=52)

## details for optimizer
opt <- list(stop = list(tol=1e-5, niter=200),
            regression = list(method="nlminb"),
            variance = list(method="nlminb"))

##########################
## models 
# A0
cntrl_A0 <- list(ar = list(f = ~ -1),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose = 1)
summary(res_A0 <- hhh4(fluBYBW,cntrl_A0))

# B0
cntrl_B0 <- list(ar = list(f = ~ 1),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_B0 <- hhh4(fluBYBW,cntrl_B0)               
 

# C0
cntrl_C0 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_C0 <- hhh4(fluBYBW,cntrl_C0)               


#A1

# weight matrix w_ji = 1/(No. neighbors of j) if j ~ i, and 0 otherwise
wji <- neighbourhood(fluBYBW)/rowSums(neighbourhood(fluBYBW))

cntrl_A1 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_A1 <- hhh4(fluBYBW,cntrl_A1)               


# B1
cntrl_B1 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_B1 <- hhh4(fluBYBW,cntrl_B1)               


# C1
cntrl_C1 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_C1 <- hhh4(fluBYBW,cntrl_C1)               


#A2
cntrl_A2 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights=wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_A2 <- hhh4(fluBYBW,cntrl_A2)               


# B2
cntrl_B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights =wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1)
res_B2 <- hhh4(fluBYBW,cntrl_B2)               

# C2
cntrl_C2 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights =wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt, verbose=1,
                 start=list(fixed=fixef(res_B0),random=c(rep(0,140),
                         ranef(res_B0)), sd.corr=c(-.5,res_B0$Sigma.orig,0)))
res_C2 <- hhh4(fluBYBW,cntrl_C2)               


# D
cntrl_D <- list(ar = list(f = ~ 1),
                ne = list(f = ~ -1 + ri(type="iid"), weights = wji),
                end = list(f =addSeason2formula(f = ~ -1 + ri(type="car") + 
                                             I((t-208)/100), S=3, period=52), 
                          offset = population(fluBYBW)),
                family = "NegBin1", optimizer = opt, verbose=1)
res_D <- hhh4(fluBYBW,cntrl_D)               

}
}

\references{
  Held, L., \enc{Höhle}{Hoehle}, M., Hofmann, M. (2005):
    A statistical framework for the analysis of multivariate infectious
    disease surveillance counts.
    \emph{Statistical Modelling}, \bold{5}, 187--199.
    
  Paul, M., Held, L. and Toschke, A. M. (2008):
    Multivariate modelling of infectious disease surveillance data.
    \emph{Statistics in Medicine}, \bold{27}, 6250--6267.    

  Paul, M. and Held, L. (2011):
    Predictive assessment of a non-linear random effects model for
    multivariate time series of infectious disease counts.
    \emph{Statistics in Medicine}, \bold{30}, 1118--1136

  Held, L. and Paul, M. (2012):
    Modeling seasonality in space-time infectious disease surveillance data.
    \emph{Biometrical Journal}, \bold{54}, 824--843
    
  Meyer, S. and Held, L. (2014):
    Power-law models for infectious disease spread.
    \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.\cr
    DOI-Link: \url{http://dx.doi.org/10.1214/14-AOAS743}

  Meyer, S., Held, L. and H\enc{ö}{oe}hle, M. (2016):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}. \emph{Journal of Statistical Software}. In press.\cr
  Preprint available at \url{http://arxiv.org/abs/1411.0416}
}

\keyword{ts}
\keyword{regression}
