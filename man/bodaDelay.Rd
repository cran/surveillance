\name{bodaDelay}
\alias{bodaDelay}
\title{Bayesian aberration detection in presence of reporting delays.}
\usage{
bodaDelay(sts, control = list(range = NULL, b = 3, w = 3,
                                          mc.munu=100, mc.y=10,
                                          pastAberrations = FALSE, 
                                          verbose = FALSE,
                                          alpha = 0.01, trend = TRUE,
                                          limit54=c(5,4), 
                                          inferenceMethod=c("asym","INLA"),
                                          noPeriods = 1, pastWeeksNotIncluded = 26,
                                          delay = TRUE))
}
\arguments{
\item{sts}{sts-object to be analysed. Needs to have a reporting triangle.}

\item{control}{list with control arguments}

\item{b}{How many years back in time to include when forming the base counts.}

\item{w}{Window's half-size, i.e. number of weeks to include before and after the current week in each year.}

\item{range}{Specifies the index of all timepoints which should be tested. If range is \code{NULL} all possible timepoints are used.}

\item{pastAberrations}{Boolean indicating whether to include an effect for past outbreaks
in a second fit of the model. This option only makes sense if \code{inferenceMethod} is \code{INLA}, as it is not supported by the other inference method.}

\item{verbose}{Boolean specifying whether to show extra debugging information.}

\item{alpha}{An approximate (one-sided) \eqn{(1-\alpha)\cdot 100\%}
prediction interval is calculated unlike the original method where it was a two-sided interval. The upper limit of this interval
i.e. the \eqn{(1-\alpha)\cdot 100\%} quantile serves as an upperbound.}

\item{trend}{Boolean indicating whether a trend should be included}

\item{noPeriods}{Number of levels in the factor allowing to use more baseline. If
equal to 1 no factor variable is created, the set of reference values is defined as in
Farrington et al (1996).}

\item{inferenceMethod}{Which inference method used, as defined in Salmon et al (2015). If one chooses \code{INLA} then inference is performed with INLA. If one chooses \code{asym} then the asymptotic normal approximation of the posteriori is used.}

\item{pastWeeksNotIncluded}{Number of past weeks to ignore in the calculation.}

\item{delay}{Boolean indicating whether to take reporting delays into account.}

\item{mc.munu}{Number of samples for the parameters of the negative binomial distribution
when performing Monte Carlo to calculate a threshold}

\item{mc.y}{Number of samples for observations
when performing Monte Carlo to calculate a threshold}

\item{limit54}{c(cases,period) is a vector allowing the
user to change these numbers.}

}
\description{
The function takes \code{range} values of the surveillance time
series \code{sts} and for each time point uses a Bayesian model of the negative binomial family with
log link inspired by the work of Noufaily et al. (2012) and of Manitz and \enc{H?hle}{Hoehle} (2014). It allows delay-corrected aberration detection as explained in Salmon et al. (2015). A \code{reportingTriangle} has to be provided in the \code{control} slot.
}
\examples{
data(salmAllOnset)
rangeTest <- 410:412
alpha <- 0.05
# Control slot for the proposed algorithm with D=0 correction
controlNormal <- list(range = rangeTest, b = 4, w = 3,
                      pastAberrations = TRUE, mc.munu=10, mc.y=10,
                      verbose = FALSE,
                      alpha = alpha, trend = TRUE,
                      limit54=c(0,50),
                      noPeriods = 10, pastWeeksNotIncluded = 26,
                      delay=FALSE,inferenceMethod="asym")
# Control slot for the proposed algorithm with D=10 correction
controlDelay <-  list(range = rangeTest, b = 4, w = 3,
                      pastAberrations = TRUE, mc.munu=10, mc.y=10,
                      verbose = FALSE,
                      alpha = alpha, trend = TRUE,
                      limit54=c(0,50),
                      noPeriods = 10, pastWeeksNotIncluded = 26,
                      delay=TRUE,inferenceMethod="asym")
salm.Normal <- bodaDelay(salmAllOnset, controlNormal)
salm.delay <- bodaDelay(salmAllOnset, controlDelay)
plot(salm.Normal)
plot(salm.delay)

}
\references{
A statistical algorithm for the early detection of outbreaks of infectious disease, Farrington, C.P., Andrews, N.J, Beale A.D. and Catchpole, M.A. (1996), J. R. Statist. Soc. A, 159, 547-563.
An improved algorithm for outbreak detection in multiple surveillance systems, Noufaily, A., Enki, D.G., Farrington, C.P., Garthwaite, P., Andrews, N.J., Charlett, A. (2012), Statistics in Medicine, published online.
Bayesian outbreak detection in the presence of reporting delays, Salmon, M., Schumacher, D., Stark, K., \enc{H?hle}{Hoehle}, M. (2015), in revision.
}
