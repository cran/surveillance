\encoding{latin1}
\name{algo.hmm}
\alias{algo.hmm}
\title{Hidden Markov Model (HMM) method}
\description{ 
  This function implements on-line HMM detection of outbreaks based on
  the retrospective procedure described in Le Strat and Carret (1999).
  Using the function \code{\link[msm]{msm}} (from package \pkg{msm}) a specified HMM
  is estimated, the decoding problem, i.e. the most probable state
  configuration, is found by the Viterbi algorithm and the most
  probable state of the last observation is recorded. On-line
  detection is performed by sequentially repeating this procedure.
  
  Warning: This function can be very slow - a more efficient implementation would be nice!
}
\usage{
  algo.hmm(disProgObj, control = list(range=range, Mtilde=-1, 
           noStates=2, trend=TRUE, noHarmonics=1,
           covEffectEqual=FALSE, saveHMMs = FALSE, extraMSMargs=list()))
}
\arguments{
        \item{disProgObj}{object of class disProg (including the observed and the state chain)}
        \item{control}{control object: 
            \describe{
            \item{\code{range}}{determines the desired time points
                which should be evaluated. Note that opposite to other
                surveillance methods an initial parameter estimation
                occurs in the HMM. Note that range should be high
                enough to allow for enough reference values for
                estimating the HMM}
            \item{\code{Mtilde}}{number of observations back in time
                to use for fitting the HMM (including
                the current observation). Reasonable values are a multiple of
                \code{disProgObj$freq}, the default is
                \code{Mtilde=-1}, which means to use all possible
                values - for long series this might take very long time!}
            \item{\code{noStates}}{number of hidden states in the HMM
                -- the typical choice is 2. The initial rates are set
                such that the \code{noStates}'th state is the one
                having the highest rate. In other words: this state is considered
                the outbreak state.}
            \item{\code{trend}}{Boolean stating whether a linear time trend exists, i.e. if \code{TRUE} (default) then \eqn{\beta_j \neq 0}{\beta != 0}}
            \item{\code{noHarmonics}}{number of harmonic waves to include in the linear predictor. Default is 1.}
            \item{\code{covEffectEqual}}{see details}
            \item{\code{saveHMMs}}{Boolean, if \code{TRUE} then the result of the fitted HMMs is saved. With this option the function can also be used to analyse data retrospectively. Default option is \code{FALSE}}
            \item{\code{extraMSMArgs}}{A named list with additional arguments to send to the \code{\link[msm:msm]{msm}} HMM fitting function. Note that the \code{msm} arguments \code{formula}, \code{data}, \code{qmatrix}, \code{hmodel}, \code{hcovariates} and \code{hconstraint} are automatically filled by \code{algo.hmm}, thus these should NOT be modified.}
            }
          }
}
\value{
  \code{algo.hmm} gives a list of class \code{survRes} which includes the 
  vector of alarm values for every timepoint in \code{range}. No
  \code{upperbound} can be specified and is put equal to zero.
    
  The resulting object contains a slot \code{control$hmm}, which
  contains the \code{msm} object with the fitted HMM.
}
\details{
  For each time point t the reference values values are extracted. If
  the number of requested values is larger than the number of possible
  values the latter is used. Now the following happens on these reference values:

  A \code{noState}-State Hidden Markov Model (HMM) is used based on
  the Poisson distribution with linear predictor on the log-link
  scale. I.e.  \deqn{Y_t | X_t = j \sim Po(\mu_t^j),}{Y_t|X_t = j ~
    Po(\mu_t^j),} where \deqn{\log(\mu_t^j) = \alpha_j + \beta_j\cdot
    t + \sum_{i=1}^{nH} \gamma_j^i \cos(2i\pi/freq\cdot (t-1)) +
    \delta_j^i \sin(2i\pi/freq\cdot (t-1))}{%
    log(mu_t^j) = alpha_j + beta_j t + \sum_{i=1}^{nH} gamma_j^i
    \cos(2*i*pi/freq * (t-1)) + delta_j^i sin(2*i*pi/freq * (t-1)) }
  and \eqn{nH=}\code{noHarmonics} and \eqn{freq=12,52} depending on the
  sampling frequency of the surveillance data. In the above \eqn{t-1} is
  used, because the first week is always saved as \code{t=1}, i.e. we
  want to ensure that the first observation corresponds to cos(0) and
  sin(0).

  If \code{covEffectEqual} then all covariate effects parameters are
  equal for the states, i.e. \eqn{\beta_j=\beta, \gamma_j^i=\gamma^i,
    \delta_j^i=\delta^i} for all \eqn{j=1,...,noState}.

  In case more complicated HMM models are to be fitted it is possible to
  modify the \code{msm} code used in this function. Using
  e.g. \code{AIC} one can select between different models (see the
  \pkg{msm} package for further details). 

  Using the Viterbi algorithms the most probable state configuration
  is obtained for the reference values and if the most probable
  configuration for the last reference value (i.e. time t) equals
  \code{control$noOfStates} then an alarm is given.

  Note: The HMM is re-fitted from scratch every time, sequential
  updating schemes of the HMM would increase speed considerably!  A
  major advantage of the approach is that outbreaks in the reference
  values are handled automatically.
}
\seealso{\code{\link[msm:msm]{msm}}}
\author{M. \enc{Höhle}{Hoehle}}
\examples{
#Simulate outbreak data from HMM
set.seed(123)
counts <- sim.pointSource(p = 0.98, r = 0.8, length = 3*52,
                              A = 1, alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.5)

\dontrun{
#Do surveillance using a two state HMM without trend component and
#the effect of the harmonics being the same in both states. A sliding
#window of two years is used to fit the HMM
surv <- algo.hmm(counts, control=list(range=(2*52):length(counts$observed),
                                   Mtilde=2*52,noStates=2,trend=FALSE,
                                   covEffectsEqual=TRUE,extraMSMargs=list()))
plot(surv,legend=list(x="topright"))
}

if (require("msm")) {
#Retrospective use of the function, i.e. monitor only the last time point
#but use option saveHMMs to store the output of the HMM fitting
surv <- algo.hmm(counts,control=list(range=length(counts$observed),Mtilde=-1,noStates=2,
                          trend=FALSE,covEffectsEqual=TRUE, saveHMMs=TRUE))

#Compute most probable state using the viterbi algorithm - 1 is "normal", 2 is "outbreak".
viterbi.msm(surv$control$hmm[[1]])$fitted

#How often correct?
tab <- cbind(truth=counts$state + 1 ,
             hmm=viterbi.msm(surv$control$hmm[[1]])$fitted)
table(tab[,1],tab[,2])
}
}
\references{
  Y. Le Strat and F. Carrat, Monitoring Epidemiologic Surveillance Data using Hidden Markov Models (1999), Statistics in Medicine, 18, 3463--3478


  I.L. MacDonald and W. Zucchini, Hidden Markov and Other Models for Discrete-valued Time Series,  (1997), Chapman & Hall, Monographs on Statistics and applied Probability 70

}

\keyword{classif}
