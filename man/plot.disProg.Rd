\name{plot.disProg}
\alias{plot.disProg}
\title{Plot Observed Counts and Defined Outbreak States of a 
  (Multivariate) Time Series}
\description{
  Plotting a (multivariate) \code{disProg} object (soft-deprecated).
  
  As of \pkg{surveillance} 1.20.0, legacy \code{disProg} objects are
  plotted via internal \code{\link{disProg2sts}} conversion
  and \code{\link{stsplot_time}}.
}
\usage{
  \method{plot}{disProg}(x, title = "", xaxis.years=TRUE, startyear = x$start[1],
    firstweek = x$start[2], as.one=TRUE, same.scale=TRUE, ...)
}
\arguments{
  \item{x}{object of class \code{disProg}}
  \item{title}{plot title}
  \item{xaxis.years}{if \code{TRUE}, the x axis is labeled using years}
  \item{startyear,firstweek}{(legacy arguments, ignored with a warning)}
  \item{as.one}{if \code{TRUE} all individual time series are shown in 
    one plot}
  \item{same.scale}{if \code{TRUE} all plots have same scale}
  \item{\dots}{further arguments passed to \code{\link{stsplot_time}}}
}
\examples{
# Plotting of simulated data
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208,
                              A = 1, alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 5)
plot(disProgObj)
title <- "Infection Counts and Defined Outbreaks for Simulated Data"
plot(disProgObj, title = title)
plot(disProgObj, title = title, xaxis.years = FALSE)

# Plotting of measles data
data(measles.weser)
# one plot
plot(measles.weser, title = "measles cases in the district Weser-Ems")
# plot cases for each "Kreis" 
plot(measles.weser, as.one = FALSE)
plot(measles.weser, as.one = FALSE, same.scale = FALSE)
}
\keyword{hplot}
\keyword{internal}
