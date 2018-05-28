\name{plot.disProg}
\alias{plot.disProg}
\alias{plot.disProg.one}
\encoding{latin1}
\title{Plot Generation of the Observed and the Defined Outbreak States of a 
  (Multivariate) Time Series}
\description{
  Plotting a (multivariate) \code{disProg} object. The internal
  function \code{plot.disProg.one} is used as a helper function
  to plot a univariate time series.
}
\usage{
  \method{plot}{disProg}(x, title = "", xaxis.years=TRUE, startyear = x$start[1],
    firstweek = x$start[2], as.one=TRUE, same.scale=TRUE, ...)

  plot.disProg.one(x, title = "", xaxis.years=TRUE, quarters=TRUE, 
    startyear =  x$start[1], firstweek = x$start[2], ylim=NULL, xlab="time",
    ylab="No. infected",type="hh",lty=c(1,1),col=c(1,1), 
    outbreak.symbol = list(pch=3, col=3), legend.opts=list(x="top", 
    legend=c("Infected", "Outbreak"), lty=NULL,pch=NULL,col=NULL), ...)
}
\arguments{
  \item{x}{object of class \code{disProg}}
  \item{title}{plot title}
  \item{xaxis.years}{if \code{TRUE}, the x axis is labeled using years}
  \item{quarters}{add quarters to the plot}
  \item{startyear}{year to begin the axis labeling (the year where
    the oldest data come from). This arguments will be obsolete in \code{sts}.}
  \item{firstweek}{number of the first week of January in the first year 
    (just for axis labeling grounds)}
  \item{as.one}{if \code{TRUE} all individual time series are shown in 
    one plot}
  \item{same.scale}{if \code{TRUE} all plots have same scale}
  \item{ylim}{range of y axis}
  \item{xlab}{label of the x-axis}
  \item{ylab}{label of the y-axis}
  \item{type}{line type of the observed counts (should be \code{hh})}
  \item{lty}{line type of the observed counts}
  \item{col}{color of the observed count lines}
  \item{outbreak.symbol}{list with entries \code{pch} and \code{col} specifying the plot symbol}
  \item{legend.opts}{a list containing the entries to be sent to the \code{\link{legend}} function. If no legend is requested use \code{legend.opts=NULL}. Otherwise, the following arguments are default
    \describe{
      \item{\code{x}}{\code{top}}
      \item{\code{legend}}{The names infected and outbreak}
      \item{\code{lty}}{If \code{NULL} the \code{lty} argument will be used}
      \item{\code{pch}}{If \code{NULL} the \code{pch} argument is used}
      \item{\code{col}}{If \code{NULL} the \code{col} argument is
        used}  
    }
    
    An further arguments to the \code{legend} function are
    just provided as additional elements of this list,
    e.g. \code{horiz=TRUE}.  
  }
  \item{\dots}{arguments passed to \code{plot.disProg.one}.
    From there, further arguments are passed to \code{\link{matplot}}.}
}
\value{
  a plot showing the number of infected and the defined alarm status
  for a time series created by simulation or given in data 
  either in one single plot or in several plots for each individual
  time series.
}
\author{M. \enc{Höhle}{Hoehle} with contributions by A. Riebler and C. Lang}
\examples{
# Plotting of simulated data
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208,
                              A = 1, alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 5)
# plot the simulated disease with the defined outbreaks
plot(disProgObj)
title <- "Number of Infected and Defined Outbreak Positions for Simulated Data"
plot(disProgObj, title = title)
plot(disProgObj, title = title, xaxis.years=TRUE,
     startyear = 1999, firstweek = 13)
plot(disProgObj, title = title, xaxis.years=TRUE,
     startyear = 1999, firstweek = 14)

# Plotting of measles data
data(measles.weser)
# one plot
plot(measles.weser, title = "measles cases in the district Weser-Ems",
     xaxis.years=TRUE, startyear= 2001, firstweek=1)
# plot cases for each "Kreis" 
plot(measles.weser, same.scale=TRUE, as.one=FALSE)
}
\keyword{hplot}
