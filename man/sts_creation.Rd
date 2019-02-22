\name{sts_creation}
\alias{sts_creation}
\title{Simulate Count Time Series with Outbreaks}
\usage{
sts_creation(theta, beta, gamma1, gamma2, m, overdispersion, dates,
  sizesOutbreak, datesOutbreak, delayMax, alpha, densityDelay)
}
\arguments{
\item{theta}{baseline frequency of reports}
\item{beta}{time trend}
\item{gamma1}{seasonality}
\item{gamma2}{seasonality}
\item{m}{seasonality}
\item{overdispersion}{\code{size} parameter of \code{\link{rnbinom}} for
  the parameterization with mean and dispersion}
\item{dates}{dates of the time series}
\item{sizesOutbreak}{sizes of all the outbreaks (vector)}
\item{datesOutbreak}{dates of all the outbreaks (vector)}
\item{delayMax}{maximal delay in time units}
\item{alpha}{alpha for getting the (1-alpha) quantile of the negative
  binomial distribution at each timepoint}
\item{densityDelay}{density distribution for the delay}
}
\description{
Function for simulating a time series and creating an
\code{\linkS4class{sts}} object.
As the counts are generated using a negative binomial distribution
one also gets the (1-alpha) quantile for each timepoint (can be interpreted
as an in-control upperbound for in-control values).
The baseline and outbreaks are created as in Noufaily et al. (2012).
}
\examples{
set.seed(12345)
# Time series parameters
scenario4 <- c(1.6,0,0.4,0.5,2)
theta <- 1.6
beta <- 0
gamma1 <-0.4
gamma2 <- 0.5
overdispersion <- 1
m <- 1
# Dates
firstDate <- "2006-01-01"
lengthT=350
dates <- as.Date(firstDate) + 7 * 0:(lengthT - 1)
# Maximal delay in weeks
D=10
# Dates and sizes of the outbreaks
datesOutbreak <- as.Date(c("2008-03-30","2011-09-25"))
sizesOutbreak <- c(2,5)
# Delay distribution
data("salmAllOnset")
in2011 <- which(isoWeekYear(epoch(salmAllOnset))$ISOYear == 2011)
rT2011 <- salmAllOnset@control$reportingTriangle$n[in2011,]
densityDelay <- apply(rT2011,2,sum, na.rm=TRUE)/sum(rT2011, na.rm=TRUE)
# alpha for the upperbound
alpha <- 0.05
# Create the sts with the full time series
stsSim <- sts_creation(theta=theta,beta=beta,gamma1=gamma1,gamma2=gamma2,m=m,
                       overdispersion=overdispersion,
                       dates=dates,
                       sizesOutbreak=sizesOutbreak,datesOutbreak=datesOutbreak,
                       delayMax=D,densityDelay=densityDelay,
                       alpha=alpha)
plot(stsSim)
}
\references{
  Noufaily, A., Enki, D.G., Farrington, C.P., Garthwaite, P., Andrews,
  N.J., Charlett, A. (2012): An improved algorithm for outbreak
  detection in multiple surveillance systems. Statistics in Medicine,
  32 (7), 1206-1222.
}
