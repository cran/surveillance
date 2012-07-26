% hhh4
\name{hhh4}
\alias{hhh4}
% functions without documentation atm
\alias{plot.ah4}
\alias{predict.ah4}
\alias{permutationTest}
\alias{pit}

\encoding{latin1}
\title{Random effects HHH model fit as described in Paul and Held (2011)} 
\description{
  Fits a Poisson or negative binomial model with conditional mean
  \deqn{\mu_{it} = \lambda_{it} y_{i,t-1} + 
                   \phi_{it} \sum_{j\neq i} w_{ji} y_{j,t-1} +
                   e_{it} \nu_{it}  }{%
        \mu_it = \lambda_it y_i,t-1 + 
                 \phi_it sum_(j != i) w_ji y_j,t-1 +
                   e_it \nu_it  }
  containing epidemic and endemic components to a multivariate time 
  series of counts. Here, the three unknown quantities 
  \itemize{
  \item \eqn{\lambda_{it}}{\lambda_it} in the autoregressive (\code{ar}) component, 
  \item \eqn{\phi_{it}}{\phi_it} in the neighbor-driven (\code{ne}) component, and
  \item \eqn{\nu_{it}}{\nu_it} in the endemic (\code{end}) component,
  }
  are decomposed additively on the log scale,
  \eqn{w_{ji}}{w_ji} are known weights, 
  and \eqn{e_{it}}{e_it} is a (multiplicative) offset.
  The linear predictors may contain random effects. 
}
\usage{
hhh4(stsObj, 
     control = list(
               ar = list(f = ~ -1,        
                         lag = 1,         
                         weights = NULL,  
                         initial = NULL   
                         ),
               ne = list(f = ~ -1,        
                         lag = 1,         
                         weights = NULL,  
                         initial = NULL   
                         ),
               end = list(f = ~ 1,         
                          offset = NULL,   
                          initial = NULL   
                          ),
               family = c("Poisson","NegBin1","NegBinM")[1],
               subset = 2:nrow(stsObj),             
               optimizer = list(tech = "nlminb",    
                                stop.tol = 1e-5,
                                stop.niter = 100),
               verbose = FALSE,
               start = list(fixed=NULL, random=NULL, sd.corr=NULL),
               data = data.frame(t = epoch(stsObj)-1),
               keep.terms = FALSE
               )
   )
}
\arguments{
  \item{stsObj}{object of class \code{"\linkS4class{sts}"} containing the multivariate
    count data time series}
  \item{control}{control object, which is a list containing several
    components:
    \describe{
      \item{\code{ar}}{ Model for the autoregressive component given as
          list with the following components: 
          \describe{
	  \item{f = ~ -1}{a formula specifying \eqn{\log(\lambda_{it})}{log(\lambda_it)}}
	  \item{lag = 1}{autoregression on \eqn{y_{i,t-lag}}{y_i,t-lag}
	    (currently only \code{lag = 1} is possible)}
	  \item{weights = NULL}{optional weights, only used if model is
	    a contact matrix (currently not implemented)}
	  \item{initial = NULL}{vector with initial values for
	    parameters if \code{f = ~1} (not really
	    used ATM)}
	}
      }
      \item{\code{ne}}{Model for the neighbor-driven component given as
       list with the following components:
	\describe{
	  \item{f = ~ -1}{a formula specifying \eqn{\log(\phi_{it})}{log(\phi_it)}}
	  \item{lag = 1}{autoregression on \eqn{y_{j,t-lag}}{y_j,t-lag}  
	    (currently only \code{lag = 1} is possible)}
	  \item{weights = NULL}{weights \eqn{w_{ji}}{w_ji}, if 
	    \code{NULL} the neighbourhood matrix of \code{stsObj} is used}
	  \item{initial = NULL}{vector with initial values for
	    parameter if \code{f = ~1} (not really used ATM)}
	}
      }
      \item{\code{end}}{Model for the endemic component given as list
	with the following components
	\describe{
	  \item{f = ~ 1}{a formula specifying \eqn{\log(\nu_{it})}{log(\nu_it)}}
	  \item{offset = NULL}{optional multiplicative offset \eqn{e_{it}}{e_it}}
	  \item{initial = NULL}{vector with initial values for
	    parameter if \code{f = ~1} (not really used ATM)}
	}
      }
      \item{\code{family}}{Distributional family -- either Poisson,
	or the negative Binomial distribution with one common
	overdispersion parameter for all units (\code{NegBin1}) or a
	overdispersion parameter for each unit (\code{NegBinM}).
      }
      \item{\code{subset}}{Typically \code{2:nrow(obs)} if model contains
	autoregression}
      \item{\code{optimizer}}{ Details for used optimizer (\code{nlminb} is 
        default, the methods specified in the \code{method} argument of \code{\link{optim}} may also be used)}
      \item{\code{verbose}}{Logical if additional information is to be printed
	during the computations}
      \item{\code{start}}{List with initials, overrides any initial
	values in formulas (this is currently the only way to specify
	initial values)}
      \item{\code{data}}{ data.frame or named list with covariates that
	are specified in the formulas for the 3 components}
      \item{\code{keep.terms}}{ Logical if the terms object used in the fit
    is to be returned }
    }
  }
}

\value{Returns an object of class \code{ah4} with elements
  \item{coefficients}{named vector with estimated (regression) parameters of the model}
  \item{se}{estimated standard errors (for regression parameters)}
  \item{cov}{covariance matrix (for regression parameters)}
  \item{Sigma}{estimated variance components for random effects}
  \item{Sigma.orig}{estimated variance components on internal scale used for optimization  }
  \item{Sigma.cov}{ inverse of (minus Hessian matrix) for variance components }
  \item{call}{ the matched call }
  \item{dim}{ vector with number of fixed and random effects in the model }
  \item{loglikelihood}{(penalized) loglikelihood evaluated at the MLE}
  \item{margll}{ (approximate) log marginal likelihood should the model contain random effects  }
  \item{convergence}{logical. Did optimizer converge?}
  \item{fitted.values}{fitted mean values \eqn{\mu_{i,t}}{\mu_it}}
  \item{control}{control object of the fit}
  \item{terms}{ the terms object used in the fit if \code{keep.terms = TRUE}
    and \code{NULL} otherwise}
  \item{stsObj}{ the supplied \code{stsObj} }
  \item{lag}{ specified lag of the autoregression, ATM always \eqn{= 1} }
  \item{nObs}{number of observations used for fitting the model}
  \item{nTime}{ number of time points used for fitting the model }
  \item{nUnit}{ number of units (e.g. areas) used for fitting the model }
}

\details{
  For further details see \code{vignette("hhh4")} and the
  references.
}

\note{
  For the time being this function performs only modelling of
  the multivariate time series. No surveillance/monitoring functionality
  is available.
}
\seealso{\code{\link{algo.hhh}}, 
\code{\link{fe}}, \code{\link{ri}}
}
\author{M. Paul and L. Held}
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
f.end2 <- addSeason2formula(f = ~ -1 + fe(1, which=rep(TRUE, ncol(measles))) + t,
                            S = 1, period = 52)
model2 <- list(ar = list(f = ~ 1), 
               end = list(f = f.end2, offset = population(measles)),
               family = "NegBin1")
# run model
summary(hhh4(measles, control=model2), idx2Exp=1, amplitudeShift=TRUE)

# include autoregressive parameter phi for adjacent "Kreise"
# no linear trend in endemic component
f.end3 <- addSeason2formula(f = ~ -1 + fe(1, which=rep(TRUE, ncol(measles))), 
                            S = 1, period = 52)
model3 <- list(ar = list(f = ~ 1),
               ne = list(f = ~1),
               end = list(f = f.end3, offset= population(measles)),
               family = "NegBin1")
# run model
summary(hhh4(measles, control=model3), idx2Exp=1:2, amplitudeShift=TRUE)

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
opt <- list(tech = "nlminb", stop.tol = 1e-5, stop.niter = 200)

##########################
## models 
# A0
cntrl_A0 <- list(ar = list(f = ~ -1),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose = 1, data = data.frame(t = epoch(fluBYBW)-1))
summary(res_A0 <- hhh4(fluBYBW,cntrl_A0))

# B0
cntrl_B0 <- list(ar = list(f = ~ 1),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_B0 <- hhh4(fluBYBW,cntrl_B0)               
 

# C0
cntrl_C0 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_C0 <- hhh4(fluBYBW,cntrl_C0)               


#A1

# weight matrix w_ji = 1/(No. neighbors of j) if j ~ i, and 0 otherwise
wji <- neighbourhood(fluBYBW)/rowSums(neighbourhood(fluBYBW))

cntrl_A1 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_A1 <- hhh4(fluBYBW,cntrl_A1)               


# B1
cntrl_B1 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_B1 <- hhh4(fluBYBW,cntrl_B1)               


# C1
cntrl_C1 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 ne = list(f = ~ 1, weights = wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_C1 <- hhh4(fluBYBW,cntrl_C1)               


#A2
cntrl_A2 <- list(ar = list(f = ~ -1),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights=wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_A2 <- hhh4(fluBYBW,cntrl_A2)               


# B2
cntrl_B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights =wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_B2 <- hhh4(fluBYBW,cntrl_B2)               

# C2
cntrl_C2 <- list(ar = list(f = ~ -1 + ri(type="iid", corr="all")),
                 ne = list(f = ~ -1 + ri(type="iid",corr="all"), weights =wji),
                 end = list(f =f.end, offset = population(fluBYBW)),
                 family = "NegBin1", optimizer = opt,
                 verbose=1, data=data.frame(t=epoch(fluBYBW)-1),
                 start=list(fixed=fixef(res_B0),random=c(rep(0,140),
                         ranef(res_B0)), sd.corr=c(-.5,res_B0$Sigma.orig,0)))
res_C2 <- hhh4(fluBYBW,cntrl_C2)               


# D
cntrl_D <- list(ar = list(f = ~ 1),
                ne = list(f = ~ -1 + ri(type="iid"), weights = wji),
                end = list(f =addSeason2formula(f = ~ -1 + ri(type="car") + 
                                             I((t-208)/100), S=3, period=52), 
                          offset = population(fluBYBW)),
                family = "NegBin1", optimizer = opt,
                verbose=1, data=data.frame(t=epoch(fluBYBW)-1))
res_D <- hhh4(fluBYBW,cntrl_D)               

}
}
\keyword{ts}
\keyword{regression}
\references{
  Held, L., \enc{Höhle}{Hoehle}, M., Hofmann, M. (2005) A statistical framework 
    for the analysis of multivariate infectious disease surveillance counts, 
    Statistical Modelling, \bold{5}, 187--199.
    
  Paul, M., Held, L. and Toschke, A. M. (2008) Multivariate modelling of 
    infectious disease surveillance data, Statistics in Medicine, \bold{27}, 
    6250--6267.    

  Paul, M. and Held, L. (2011) Predictive assessment of a non-linear
    random  effects model for multivariate time series of infectious
    disease counts. Statistics in Medicine, \bold{30}, 1118--1136
}




