\encoding{latin1}
\name{hhh4_formula}
\alias{fe}
\alias{ri}

\title{
Specify Formulae in a Random Effects HHH Model
}

\description{
The special functions \code{fe} and \code{ri} are used to specify 
(unit-specific) effects of covariates and a random intercept term, 
respectively, in formulae used in the function \code{\link{hhh4}}.
}

\usage{
  fe(x, which = NULL, initial = NULL)

  ri(type = c("iid","car")[1], corr = c("none", "all")[1],
     initial.var = NULL, initial.re = NULL)
}

\arguments{
  \item{x}{name of the covariate, use 1 for an intercept}
  \item{which}{vector of logicals indicating which unit
   should get an unit-specific parameter. If \code{NULL}, 
   the effect of the covariate is assumbed to be the same
   for all units.}
  \item{initial}{initial values (on internal scale!) 
    for the fixed effects used for optimization.
    Not really used ATM.}
  \item{type}{random intercepts either follow an IID or a CAR
   model}
  \item{corr}{logical switch indicating whether random effects
   in different components (such as \code{ar} and \code{end})
   should be correlated or not. }
  \item{initial.var}{initial values (on internal scale!) 
    for the variance components used for optimization. Not really used ATM.}
  \item{initial.re}{initial values (on internal scale!) for the random effects 
    used for optimization. Not really used ATM.}
}


\seealso{
\code{\link{addSeason2formula}}, 
usage of formulae in the vignette and in examples of \code{\link{hhh4}}.
}

\note{
This function should only be used in formula objects for \code{\link{hhh4}}, 
and is not intended for direct calling.

All covariates specified in \code{fe} must be available in the
element \code{data} of the \code{control} argument of \code{\link{hhh4}}.

If unit-specific or random intercepts are specified, an overall
intercept must be excluded with \code{-1}.
}

\keyword{regression}


\examples{
\dontrun{
# some calls of the fitting function 'hhh4': 
# see vignette("hhh4") for further details

## Ex: univariate time series of meningococcal infections in Germany
# Negative binomial model with
# endemic component: Intercept + S = 1 sine/cosine pair
# autoregressive component: Intercept

f.S1 <- addSeason2formula(f = ~ 1, S = 1, period = 52)
hhh4(meningo, control = list(ar = list(f = ~ 1),
                             end = list(f = f.S1), 
                             family = "NegBin1")) 


## Ex: disease-specific intercept in influenza/meningococcal time series
# Negative binomial model with
# autoregressive component: disease-specific intercepts
# neighbour-driven component: only intercept for flu -> men
# endemic component: S=3 and S=1 sine/cosine pairs for flu and men, respectively

f.end <- addSeason2formula(f = ~ -1 + fe(1, which = c(TRUE, TRUE)), 
                           S = c(3, 1),   
                           period = 52)
m <- list(ar = list(f = ~ -1 + fe(1, which = c(TRUE, TRUE))), 
          ne = list(f = ~ -1 + fe(1, which = c(FALSE, TRUE))), 
          end = list(f = f.end),
          family = "NegBinM"
          )
hhh4(fluMen, control = m)


## Ex: (correlated) random intercepts for influenza in Southern Germany
# Negative binomial model with
# autoregressive component: Intercept
# neighbour-driven component: random intercepts
# endemic component: random intercepts + trend + S = 3 sine/cosine pairs

f.end <- addSeason2formula(f = ~ -1 + ri(type = "iid", corr="all") +
                                 I((t-208)/100), S = 3, period = 52)

model.B2 <- list(ar = list(f = ~ 1),
                 ne = list(f = ~ -1 + ri(type = "iid", corr="all"), 
                           weights = wji),
                 end = list(f = f.end, offset = population(flu)),
                 family = "NegBin1")
hhh4(flu, model.B2)


## Ex: measles in Germany
# Poisson model with
# autoregressive component: Intercept + vaccination coverage info
# endemic component: Intercept + S = 1 sine/cosine pair

f.end <- addSeason2formula(f = ~ 1, S = 1, period = 26)
model.A0 <- list(ar = list(f = ~ 1 + logVac0),
                 end = list(f = f.end, offset = population(measles2w)),
                 data = list(t = epoch(measles2w), logVac0 = log(vac0)))
hhh4(measles2w, model.A0)    
}
}


