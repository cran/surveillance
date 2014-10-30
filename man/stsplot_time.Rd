\encoding{latin1}
\name{stsplot_time}

\alias{stsplot_time}
\alias{stsplot_time1}
\alias{stsplot_alarm}
% helper functions for time axis formatting
\alias{atChange}
\alias{at2ndChange}
\alias{atMedian}

\title{
  Time-Series Plots for \code{"sts"} Objects
}

\description{
  These are the \code{plot} variants of \code{type=observed~time|unit},
  \code{type=observed~time}, and \code{type=alarm~time}
  for \code{"\linkS4class{sts}"} objects (see the central \code{"sts"}
  \code{\link[=plot,sts,missing-method]{plot}}-method for
  an overview of plot types).
}

\usage{
stsplot_time(x, units, method=x@control$name, disease=x@control$data,
             as.one=FALSE, same.scale=TRUE, par.list=list(), ...)

stsplot_time1(x, k=1, ylim=NULL,
              axes=TRUE, xaxis.tickFreq=list("\%Q"=atChange),
              xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="\%G\n\n\%OQ",
              epochsAsDate=x@epochAsDate,
              xlab="time", ylab="No. infected", main=NULL,
              type="s", lty=c(1,1,2), col=c(NA,1,4), lwd=c(1,1,1),
              outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1),
              alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
              legend.opts=list(x="top",legend=NULL,lty=NULL,pch=NULL,col=NULL),
              dx.upperbound=0L, hookFunc=function(){},
              .hookFuncInheritance=function() {}, ...)

stsplot_alarm(x, lvl=rep(1,nrow(x)), ylim=NULL,
              xaxis.tickFreq=list("\%Q"=atChange),
              xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="\%G\n\n\%OQ",
              epochsAsDate=x@epochAsDate, xlab="time", main=NULL,
              type="hhs", lty=c(1,1,2), col=c(1,1,4),
              outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1),
              alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
              cex=1, cex.yaxis=1, ...)
}

\arguments{
  \item{x}{an object of class \code{"\linkS4class{sts}"}.}
  \item{units}{optional integer vector to select the units (=columns of
    \code{observed(x)}) to plot. \code{stsplot_time1} is called for
    \code{k in units}. The default is to plot all time series.}
  \item{method}{name of the surveillance method to be used in the main
    title. Currently ignored.}
  \item{disease}{name of the disease to be used in the main title.
    Currently ignored.}
  \item{as.one}{logical indicating if all time series should be plotted
    within the same frame. This is currently not implemented
    for the \code{"sts"} class, but see \code{\link{plot.disProg}} and
    \code{\link{sts2disProg}}.}
  \item{same.scale}{logical indicating if all time series should be
    plotted with the same \code{ylim}. Default is to do so. Only
    relevant for multivariate plots (\code{ncol(x) > 1}).}
  \item{par.list}{a list of arguments delivered to a call of
    \code{\link{par}} before each plot.}

  \item{k}{the unit to plot, i.e., an element of \code{1:ncol(x)}.}
  \item{ylim}{the y limits of the plot(s). Ignored if
    \code{same.scale=FALSE}.}
  \item{axes}{a logical value indicating whether both axes should be drawn
    on the plot.}
  \item{xaxis.tickFreq,xaxis.labelFreq,xaxis.labelFormat}{see Details.}
  \item{epochsAsDate}{Boolean indicating whether to treat the epochs as
    Date objects. Default: Value of the corresponding slot in \code{x}.}
  \item{xlab}{a title for the x axis. See \code{plot.default}.}
  \item{ylab}{a title for the y axis. See \code{plot.default}.}
  \item{main}{an overall title for the plot: see 'title'.}
  \item{type}{type of plot to do.}
  \item{lty}{vector of length 3 specifying the line type for the three
    lines in the plot -- see \code{col} argument.}
  \item{col}{Vector of length 3 specifying the color to use in the
    plot. The first color is the fill color of the polygons for the
    counts bars (\code{NA} for unfilled), the 2nd element denotes their
    border color, the 3rd element is the color of the \code{upperbound}
    plotting.} 
  \item{lwd}{Vector of length 3 specifying the line width of the three
    elements to plot. See also the \code{col} argment.}
  \item{alarm.symbol}{a list with entries \code{pch}, \code{col},
    \code{cex} and \code{lwd} specifying the appearance of the outbreak
    symbol in the plot.}
  \item{outbreak.symbol}{a list with entries \code{pch}, \code{col},
    \code{cex} and \code{lwd} specifying the appearance of the outbreak
    symbol in the plot.}
  \item{legend.opts}{a list containing the entries to be sent to the
  \code{\link{legend}} function. If no legend is requested use
  \code{legend.opts=NULL}. Otherwise, the following arguments are
  default 
    \describe{
      \item{\code{x}}{\code{top}}
      \item{\code{legend}}{The names infected and outbreak.}
      \item{\code{lty}}{If \code{NULL} the \code{lty} argument will be used}
      \item{\code{pch}}{If \code{NULL} the \code{pch} argument is used}
      \item{\code{col}}{If \code{NULL} the \code{col} argument is
        used}  
    }
    
    Any further arguments to the \code{legend} function are
    just provided as additional elements of this list,
    e.g. \code{horiz=TRUE}.
  }
  \item{dx.upperbound}{horizontal change in the plotting of the
    upperbound line. Sometimes it can be convenient to offset this line
    a little for better visiability.}
  \item{lvl}{A vector of length \code{ncol(x)}, which is used to
    specify the hierarchy level for each time series in the sts object
    for alarm plots.} 
  \item{cex}{A numerical value giving the amount by which plotting text
    and symbols should be magnified relative to the default. See
    \code{?par} for details.}
  \item{cex.yaxis}{The magnification to be used for y-axis annotation
    relative to the current setting of 'cex'.}
  \item{hookFunc}{a function that is called after all the basic plotting
    has be done, i.e., it is not possible to control formatting with
    this function. See Examples.}
  \item{.hookFuncInheritance}{a function which is altered by sub-classes
    plot method. Do not alter this function manually.}
  \item{...}{further arguments for the function \code{matplot}. If
    e.g. \code{xlab} or \code{main} are provided they overwrite the
    default values.} 
}

\details{
  The time series plot relies on the work-horse \code{stsplot_time1}.
  Its arguments are (almost) similiar to \code{\link{plot.survRes}}.

  In case the epochs of the \code{sts} object are \code{Date} objects it
  is possible to obtain very flexible formatting of the x-axis and its
  annotations using the \code{xaxis.tickFreq}, \code{xaxis.labelFreq}
  and \code{xaxis.labelFormat}. The first two are named lists containing
  pairs with the \emph{name} being a \code{\link{strftime}} single
  conversion specification and the second part is a function which based
  on this conversion returns a subset of the rows in the \code{sts}
  objects. The subsetting function has the following header:
  \code{function(x,xm1)}, where \code{x} is a vector containing
  the result of applying the conversion in \code{name} to the epochs of
  the \code{sts} object and \code{xm1} is the scalar result when
  applying the conversion to the natural element just before the first
  epoch.  Three predefined subsetting functions exist:
  \code{atChange}, \code{at2ndChange}  and \code{atMedian}, which 
  are used to make a tick at each (each 2nd for \code{at2ndChange})
  change and at the median index computed on all having the same value,
  respectively: 
  \preformatted{
    atChange <- function(x,xm1) which(diff(c(xm1,x)) != 0)
    at2ndChange <- function(x,xm1) which(diff(c(xm1,x) \%/\% 2) != 0)
    atMedian <- function(x,xm1) tapply(seq_along(x), INDEX=x, quantile, prob=0.5, type=3)
  }
  By defining own functions here, one can obtain an arbitrary degree of
  flexibility. 
  Finally, \code{xaxis.labelFormat} is a \code{\link{strftime}}
  compatible formatting string., e.g. the default value is
  \code{"\%G\\n\\n\%OQ"}, which means ISO year and quarter (in roman
  letters) stacked on top of each other.
  
}

\value{
  \code{NULL} (invisibly).
  The functions are called for their side-effects.
}

\author{
  Michael H\enc{ö}{oe}hle and Sebastian Meyer
}

\seealso{
  the central \code{\link{stsplot}}-documentation for an overview of
  plot types.
}

\examples{
data("ha.sts")
print(ha.sts)

plot(ha.sts, type=observed ~ time | unit)  # default multivariate type
plot(ha.sts, type=observed ~ time)         # aggregated over all districts

## Hook function example
hookFunc <- function() grid(NA,NULL,lwd=1)
plot(ha.sts, hookFunc=hookFunc)

## Use ISO8601 date formatting (see ?strptime) and no legend
data("salmNewport")
plot(aggregate(salmNewport,by="unit"), xlab="Time (weeks)",
     xaxis.tickFreq=list("\%m"=atChange,"\%G"=atChange),
     xaxis.labelFreq=list("\%G"=atMedian),xaxis.labelFormat="\%G",
     legend.opts=NULL)

## Formatting now also works for daily data (illustrate by artifical
## outbreak converted to sts object by linelist2sts)
set.seed(123)
exposureTimes <-  as.Date("2014-03-12") + sample(x=0:25,size=99,replace=TRUE)
sts <- linelist2sts(data.frame(exposure=exposureTimes),
                               dateCol="exposure",aggregate.by="1 day")
## Plot it with larger ticks for days than usual
surveillance.options("stsTickFactors"=c("\%d"=1, "\%W"=0.33,
                "\%V"=0.33, "\%m"=1.75, "\%Q"=1.25, "\%Y"=1.5, "\%G"=1.5))
plot(sts,xaxis.tickFreq=list("\%d"=atChange,"\%m"=atChange),
     xaxis.labelFreq=list("\%d"=at2ndChange),xaxis.labelFormat="\%d-\%b",
     legend.opts=NULL, xlab="Time (days)")
}

\keyword{hplot}
\keyword{ts}
