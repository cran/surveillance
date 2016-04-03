\encoding{latin1}
\name{algo.twins}
\alias{algo.twins}

\title{Model fit based on a two-component epidemic model} 
\description{
  Fits a negative binomial model (as described in Held et al. (2006)
  to an univariate time series of counts.  
}

\usage{
algo.twins(disProgObj, control=list(burnin=1000, filter=10,
   sampleSize=2500, noOfHarmonics=1, alpha_xi=10, beta_xi=10,
   psiRWSigma=0.25,alpha_psi=1, beta_psi=0.1, nu_trend=FALSE,
   logFile="twins.log"))
}

\arguments{
  \item{disProgObj}{object of class \code{disProg}}
  \item{control}{control object:
        \describe{
        \item{\code{burnin}}{Number of burn in samples.}
        \item{\code{filter}}{Thinning parameter. If \code{filter = 10} every 10th sample is after the burn in is returned.}
        \item{\code{sampleSize}}{Number of returned samples. Total
	  number of samples = \code{burnin}+\code{filter}*\code{sampleSize}}
	\item{\code{noOfHarmonics}}{Number of harmonics to use in the
	  modelling, i.e. \eqn{L}{L} in (2.2) of Held et al (2006).}
        \item{\code{alpha_xi}}{Parameter \eqn{\alpha_{\xi}}{\alpha_\xi} of the hyperprior of the epidemic parameter \eqn{\lambda}{\lambda}}
        \item{\code{beta_xi}}{Parameter \eqn{\beta_{\xi}}{\beta_\xi} of the hyperprior of the epidemic parameter \eqn{\lambda}{\lambda}}
        \item{\code{psiRWSigma}}{Starting value for the tuning of the variance of the random walk proposal for the overdispersion parameter \eqn{\psi}{\psi}.}
        \item{\code{alpha_psi}}{Parameter \eqn{\alpha_{\psi}}{\alpha_\psi} of the prior of the overdispersion parameter \eqn{\psi}{\psi}}
        \item{\code{beta_psi}}{Parameter \eqn{\beta_{\psi}}{\beta_\psi}
	  of the prior of the overdispersion parameter \eqn{\psi}{\psi}}
	\item{\code{nu_trend}}{Adjust for a linear trend in the endemic
	  part? (default: \code{FALSE})}
        \item{\code{logFile}}{Base file name for the output files. The function writes three output files in the current working directory \code{getwd()}. If \code{logfile = "twins.log"} the results are stored in the three files \file{twins.log}, \file{twins.log2} and \file{twins.log.acc}.\cr
        \file{twins.log} contains the returned samples of the parameters \eqn{\psi}{\psi}, \eqn{\gamma_{0}}{\gamma_0}, \eqn{\gamma_{1}}{\gamma_1}, \eqn{\gamma_{2}}{\gamma_2}, K, \eqn{\xi_{\lambda}}{\xi_\lambda} \eqn{\lambda_{1},...,\lambda{n}}{\lambda_1,...,\lambda_{n}}, the predictive distribution of the number of cases at time \eqn{n+1} and the deviance.\cr
        \file{twins.log2} contains the sample means of the variables \eqn{X_{t}, Y_{t}, \omega_{t}}{X_t, Y_t, \omega_t} and the relative frequency of a changepoint at time t for t=1,...,n and the relative frequency of a predicted changepoint at time n+1.\cr
        \file{twins.log.acc} contains the acceptance rates of \eqn{\psi}{\psi}, the changepoints and the endemic parameters \eqn{\gamma_{0}}{\gamma_0}, \eqn{\gamma_{1}}{\gamma_1}, \eqn{\gamma_{2}}{\gamma_2} in the third column and the variance of the random walk proposal for the update of the parameter \eqn{\psi}{\psi} in the second column.}
}
}
}

\details{Note that for the time being this function is not a
  surveillance algorithm, but only a modelling approach as described in
  the Held et. al (2006) paper.

  Note also that the function writes three logfiles in the current
  working directory \code{getwd()}: \file{twins.log},
  \file{twins.log.acc} and \file{twins.log2}.
  Thus you need to have write permissions in the current working
  directory.

  Finally, inspection of the C++ code using valgrind shows some memory
  leaks when running the old underlying C++ program. As we are unable to
  fix this impurity at the present time, we have instead put the example
  code in a 'dontrun' environment. The example code, however, works fine
  -- the measure is thus more aimed at reducing the number of CRAN
  problems with the package.
  
}

\value{Returns an object of class \code{atwins} with elements 
  \item{control}{specified control object}
  \item{disProgObj}{specified \code{disProg}-object}
  \item{logFile}{contains the returned samples of the parameters \eqn{\psi}{\psi}, \eqn{\gamma_{0}}{\gamma_0}, \eqn{\gamma_{1}}{\gamma_1}, \eqn{\gamma_{2}}{\gamma_2}, K, \eqn{\xi_{\lambda}}{\xi_\lambda} \eqn{\lambda_{1},...,\lambda{n}}{\lambda_1,...,\lambda_{n}}, the predictive distribution and the deviance.}
  \item{logFile2}{contains the sample means of the variables \eqn{X_{t}, Y_{t}, \omega_{t}}{X_t, Y_t, \omega_t} and the relative frequency of a changepoint at time t for t=1,...,n and the relative frequency of a predicted changepoint at time n+1.} 
}
\references{
  Held, L., Hofmann, M., \enc{Höhle}{Hoehle}, M. and Schmid V. (2006):
  A two-component model for counts of infectious diseases.
  \emph{Biostatistics}, \bold{7}, pp. 422--437.
}

\author{
  M. Hofmann and M. \enc{Höhle}{Hoehle} and
  D. \enc{Sabanés Bové}{Sabanes Bove}
}

\examples{
\dontrun{
# Load the data used in the Held et al. (2006) paper
data("hepatitisA")

# Fix seed - this is used for the MCMC samplers in twins
set.seed(123)

# Call algorithm and save result (use short chain without filtering for speed)
otwins <- algo.twins(hepatitisA,
                     control=list(burnin=500, filter=1, sampleSize=1000))

# This shows the entire output (use ask=TRUE for pause between plots)
plot(otwins, ask=FALSE)

# Direct access to MCMC output
hist(otwins$logFile$psi,xlab=expression(psi),main="")
if (require("coda")) {
    print(summary(mcmc(otwins$logFile[,c("psi","xipsi","K")])))
}
}
}

\keyword{ts}
\keyword{regression}
