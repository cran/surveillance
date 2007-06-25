\name{algo.farrington}
\alias{algo.farrington}
\title{Surveillance for a time series using the Farrington procedure.}
\description{
  The function takes \code{range} values of the time series \code{counts} and for each uses a GLM to predict the number of counts according to the procedure by Farrington et. al. This is then compared to the observed number of counts and in case an exceedance of the confidence interval calculated is seen an alarm is raised. 
}
\usage{
  algo.farrington(disProgObj, control=list(range=NULL, b=3, w=3,
  reweight=TRUE,verbose=FALSE,alpha=0.01))
}
\arguments{
\item{disProgObj}{object of class disProgObj (including the observed and the state chain)}
\item{control}{Control object
    \itemize{
    \item{range}{Specifies the index of all timepoints which should be tested. If range is \code{NULL} the maximum number of possible weeks is used.}
    \item{b}{how many years back in time to include when forming the base counts.}
    \item{w}{windows size, i.e. number of weeks to include before and after the current week}
    \item{reweight}{Boolean specifying whether to perform reweight step}
    \item{trend}{If \code{true} a trend is included and kept in case
        the conditions in the Farrington et. al. paper are met (see
        the results). If \code{false} then NO trend is fit.}
    \item{verbose}{show extra debugging information}
    \item{plot}{shows the final GLM model fit graphically (use History|Recording to see all pictures)}
    \item{alpha}{An approximate (two-sided) \eqn{(1-\alpha)\%}\ confidence interval is calculated}
    \item{limit54}{to avoid alarms in scenarios going from zero to 1-2 cases the algorithm uses the following criterion (see Section 3.8 of the Farrington paper): no alarm is sounded if fewer than \eqn{cases=5} reports were received in the past \eqn{period=4} weeks. \code{o54=c(cases,perior)} is a vector allowing the user to change these numbers}
    }}
}
\details{
  The following steps are perfomed according to the Farrington
  et. al. paper.
\enumerate{
\item fit of the initial model and initial estimation of mean and overdispersion.
\item calculation of the weights omega (correction for past outbreaks)
\item refitting of the model
\item revised estimation of overdispersion
\item rescaled model
\item omission of the trend, if it is not significant
\item repetition of the whole procedure
\item calculation of the threshold value
\item computation of exceedance score
}
}
\value{
An object of class \code{SurvRes}.
}

\examples{
#Read Salmonella Agona data
data("salmonella.agona")

#Do surveillance for the last 100 weeks.
n <- length(salmonella.agona$observed)
#Set control parameters.
control <- list(b=4,w=3,range=(n-100):n,reweight=TRUE, verbose=FALSE,alpha=0.01)
res <- algo.farrington(salmonella.agona,control=control)
#Plot the result.
plot(res,disease="Salmonella Agona",method="Farrington")
}
\author{M. Höhle}
\seealso{\code{\link{algo.farrington.fitGLM}},\code{\link{algo.farrington.threshold}}}
\keyword{}
\source{A statistical algorithm for the early detection of outbreaks of infectious disease, Farrington, C.P., Andrews, N.J, Beale A.D. and Catchpole, M.A. (1996). , J. R. Statist. Soc. A, 159, 547-563.}


