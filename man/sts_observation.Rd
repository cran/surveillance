\name{sts_observation}
\alias{sts_observation}
\title{Create an \code{sts} object with a given observation date}
\usage{
sts_observation(sts, dateObservation, cut = TRUE)
}
\arguments{
\item{sts}{sts-object we want to set at a previous state. Needs to include a reporting triangle.}
\item{dateObservation}{Date for which we want the state. Needs to be in the reporting triangle dates.}
\item{cut}{Boolean indicating whether to have 0 counts after the observation date or to simply cut the sts-object}
}
\description{
Function for creating an \code{\linkS4class{sts}} object with a given observation date.
}
\examples{
data("salmAllOnset")
salmAllOnsety2013m01d20 <- sts_observation(salmAllOnset,
dateObservation="2014-01-20",cut=FALSE)
plot(salmAllOnset)
lines(salmAllOnsety2013m01d20@observed,t="h",col="red")
}
