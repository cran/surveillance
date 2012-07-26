\name{display-methods}
\alias{plot,sts,missing-method}
%\alias{print}
%\alias{print,sts-method}
\alias{show,sts-method}
\alias{plot.sts.time}
\alias{plot.sts.time.one}
\alias{plot.sts.spacetime}
\alias{plot.sts.alarm}

\title{Display Methods for Surveillance Time-Series Objects}
\description{
  Display methods for objects of class \code{"sts"}.
}

\section{Usage}{
  plot(x,y,type,...)
  print(x,...)
}

\section{Arguments}{
\describe{
  \item{x}{an object of class \code{"sts"}}
  \item{y}{missing}
  \item{type}{a formula specifying the plot type, several options are
    possible:
    \describe{
	\item{\code{observed ~ time}}{The observations in \code{x} are
	  aggregated over units and the resulting univariate time-series
	  is plotted. The plotting is done by the function
	  \code{plot.time.sts}, which takes the same arguments as the
	  \code{\link{plot.survRes}} function.}
	\item{\code{observed ~ time | unit}}{shows \code{dim(x)} plots
	  with each showing the time-series of one observational
	  unit. The actual plotting is done by the function 
	  \code{plot.time.sts.one}}
	\item{\code{observed ~ 1 | unit}}{for each unit the counts are
	  aggregated over time and a map illustrating the counts is
	  shown. The column names of the \code{x@observed} object are
	  used to label the entries of the \code{x@map}. Regions with an
	  alarm are shaded.}
	\item{\code{observed ~ 1 | unit * time}}{an animation consisting
	  of \code{nrow(x)} frames is generated. Each frame contains the
	  number of counts per region for the current row in the
	  \code{observed} matrix. It is possible to redirect the output
	  into files, e.g. to generate an animated GIF. See the
	  examples.}
	\item{\code{alarm ~ time}}{Generates a so called alarmplot for a
	  multivariate \code{sts} object. For each time point and each
	  series it is shown whether there is an alarm. In case of
	  hierarchical surveillance the user can pass
	  an additional argument \code{lvl}, which is a vector of the
	  same length as rows in \code{x} specifying for each time series
	  its level.
	}
      }
    }
    \item{\dots}{further arguments passed to or from other methods:
      in case of plotting these are passed to \code{\link{plot}},
      in case of printing these are passed to \code{\link{print.default}}}
  }
}
\details{
  The plotting of time-series plots relies on two internal functions
  with \code{plot.sts.time.one} being the work-horse. Its arguments are
  (almost) similiar to \code{\link{plot.survRes}}. \code{k} is the
  column to plot.

  \preformatted{plot.sts.time(x,
    type,method=x@control$name,disease=x@control$data,same.scale=TRUE,
    par.list=list(mfrow=magic.dim(nAreas),mar=par()$mar),...)}
  
  \preformatted{plot.sts.time.one <- function(x, k=1,
  domany=FALSE,ylim=NULL,xaxis.years=TRUE, xaxis.units=TRUE,
  epochsAsDate=x@epochAsDate, xlab="time", ylab="No. infected",
  main=NULL, type="s",lty=c(1,1,2),col=c(NA,1,4),lwd=c(1,1,1),
  outbreak.symbol = list(pch=3, col=3, cex=1),alarm.symbol=list(pch=24,
  col=2, cex=1),cex=1,legend.opts=list(x="top",
  legend=NULL,lty=NULL,pch=NULL,col=NULL),dx.upperbound=0.5,
  hookFunc=function() {},...)}

  Note that the \code{hookFunc} is called after all the basic plotting
  has be done, i.e. it is not possible to control formatting with this
  function.
  
  For spacetime plots the following internal function does all the work:

  \preformatted{plot.sts.spacetime(x,type,legend=NULL,opts.col=NULL,
    labels=TRUE,wait.ms=250,cex.lab=0.7,verbose=FALSE,
    dev.printer=NULL,...)}

  Alarmplot
  \preformatted{plot.sts.alarm <- function(x, lvl=rep(1,nrow(x)),
  ylim=NULL,xaxis.years=TRUE, xaxis.units=TRUE,
  epochsAsDate=x@epochAsDate, xlab="time", main=NULL,
  type="hhs",lty=c(1,1,2),col=c(1,1,4), outbreak.symbol = list(pch=3,
  col=3, cex=1),alarm.symbol=list(pch=24, col=2,
  cex=1),cex=1,cex.yaxis=1,...)}

  \code{print} is the method for printing \code{sts} objects.

}

\value{

The methods are called for their side-effects.

}

\seealso{
  \code{\link{plot.survRes}}
}

\examples{
data("ha.sts")
print(ha.sts)
plot(ha.sts, type= observed ~ time)
plot(ha.sts, type= observed ~ time | unit)
plot(ha.sts, type= observed ~ 1 | unit)
plot(ha.sts[1:20,1:2], type= observed ~ 1 | unit)
plot(aggregate(ha.sts,nfreq=13), type= observed ~ 1 | unit * time)


\dontrun{
#Configure a png device printer for the plot command
dev.printer <-
list(device=png,extension=".png",width=640,height=480,name="/tmp/berlin")

#Do the animation
plot(aggregate(ha.sts,nfreq=13), type = observed ~ 1 | unit * time,
         dev.printer=dev.printer)

#Use ImageMagick -- replace /sw/bin/convert by your path to convert
system(paste("/sw/bin/convert -delay 50 ",dev.printer$name,
           "*.png /tmp/animated.gif",sep=""))
}

}

\keyword{ts}
\keyword{hplot}
