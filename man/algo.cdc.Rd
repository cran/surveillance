\name{algo.cdc}
\alias{algo.cdcLatestTimepoint}
\alias{algo.cdc}
\encoding{latin1}

\title{The CDC Algorithm}
\description{
  Surveillance using the CDC Algorithm
}
\usage{
algo.cdcLatestTimepoint(disProgObj, timePoint = NULL,
                        control = list(b = 5, m = 1, alpha=0.025))
algo.cdc(disProgObj, control = list(range = range, b= 5, m=1, 
         alpha = 0.025))
}
\arguments{
  \item{disProgObj}{object of class disProg (including the observed and the state chain).}
  \item{timePoint}{time point which should be evaluated in \code{algo.cdcLatestTimepoint}. The
    default is to use the latest timepoint.}
  \item{control}{control object: \code{range} determines the desired
    timepoints which should be evaluated, \code{b} describes the number of years to go
    back for the reference values, \code{m} is the half window width for the reference
    values around the appropriate timepoint (see details).
    The standard definition is \code{b}=5 and \code{m}=1.}
}
\details{
  Using the reference values for calculating an upper limit, alarm is
  given if the actual value is bigger than a computed threshold.
  \code{algo.cdc} calls \code{algo.cdcLatestTimepoint} for the values
  specified in \code{range} and for the system specified in
  \code{control}.  The threshold is calculated from the predictive
  distribution, i.e.  \deqn{mean(x) + z_{\alpha/2} * sd(x) * \sqrt(1+1/k),}
  which corresponds to Equation 8-1 in Farrington and Andrews (2003).
  Note that an aggregation into 4-week blocks occurs in
  \code{algo.cdcLatestTimepoint} and \code{m} denotes number of 4-week
  blocks (months) to use as reference values. This function currently
  does the same for monthly data (not correct!)
}
\value{
  \code{algo.cdcLatestTimepoint} returns a list of class \code{survRes} (surveillance result), which
  includes the alarm value (alarm = 1, no alarm = 0) for recognizing an
  outbreak, the threshold value for recognizing the alarm and
  the input object of class disProg.
  
  \code{algo.cdc} gives a list of class \code{survRes} which
  includes the vector of alarm values for every timepoint in
  \code{range}, the vector of threshold values for every timepoint
  in \code{range} for the system specified by \code{b}, \code{w},
  the range and the input object of class disProg.
}
\seealso{
  \code{\link{algo.rkiLatestTimepoint}},\code{\link{algo.bayesLatestTimepoint}}
  and \code{\link{algo.bayes}} for the Bayes system.
}
\author{M. \enc{Höhle}{Hoehle}}
\examples{
# Create a test object
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 500, 
                              A = 1,alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.7)

# Test week 200 to 208 for outbreaks with a selfdefined cdc
algo.cdc(disProgObj, control = list(range = 400:500,alpha=0.025))
}

\keyword{classif}
\references{
Stroup, D., G. Williamson, J. Herndon, and J. Karon (1989). Detection
of aberrations in the occurence of notifiable diseases surveillance data.
Statistics in Medicine 8, 323-329. 

Farrington, C. and N. Andrews (2003). Monitoring the Health of
Populations, Chapter Outbreak Detection: Application to Infectious
Disease Surveillance, pp. 203-231. Oxford University Press.
}
