\docType{data}
\name{salmNewportTruncated}
\alias{salmNewportTruncated}
\title{Salmonella Newport cases in Germany 2001-2011 by data of symptoms onset}
\format{A list of 6 sts-objects}
\usage{
data("salmNewportTruncated")
}

\examples{
  data("salmNewportTruncated")
  par(mfrow=c(2,3))
  for (week in 44:49){
  plot(salmNewportTruncated[[week]],ylim=c(0,62))
}
}
\description{
  A list containing 6 time series of counts of Salmonella Newport in Germany 2001-2011 aggregated by data of disease onset as known during weeks 44 to 49 of 2011. 
}
\references{
  Bayesian outbreak detection in the presence of reporting delays, Salmon, M., Schumacher, D., Stark, K., \enc{H?hle}{Hoehle}, M. (2015), in revision.
}
\keyword{datasets}
