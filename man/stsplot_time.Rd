\encoding{latin1}
\name{stsplot_time}

\alias{stsplot_time}
\alias{stsplot_time1}
\alias{stsplot_alarm}

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
stsplot_time(x, units=NULL,
             as.one=FALSE, same.scale=TRUE, par.list=list(), ...)

stsplot_time1(x, k=1, ylim=NULL,
              axes=TRUE, xaxis.tickFreq=list("\%Q"=atChange),
              xaxis.labelFreq=xaxis.tickFreq, xaxis.labelFormat="\%G\n\n\%OQ",
              epochsAsDate=x@epochAsDate,
              xlab="time", ylab="No. infected", main=NULL,
              type="s", lty=c(1,1,2), col=c(NA,1,4), lwd=c(1,1,1),
              outbreak.symbol=list(pch=3, col=3, cex=1, lwd=1),
              alarm.symbol=list(pch=24, col=2, cex=1, lwd=1),
              legend.opts=list(),
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
  \item{units}{optional integer or character vector to select the units (=columns of
    \code{observed(x)}) to plot. The default is to plot all time series.
    If \code{as.one=FALSE}, \code{stsplot_time1} is called
    \code{for (k in units)} with \code{mfrow} splitting (see \code{par.list}).
    Note that if there are too many \code{units}, the default \code{mfrow}
    setting might lead to the error \dQuote{figure margins too large}
    (meaning that the units do not fit onto a single page).}
  \item{as.one}{logical indicating if all time series should be plotted
    in a single frame (using \code{\link{matplot}}).}
  \item{same.scale}{logical indicating if all time series should be
    plotted with the same \code{ylim}. Default is to do so. Only
    relevant for multivariate plots (\code{ncol(x) > 1}).}
  \item{par.list}{a list of arguments delivered to a call of
    \code{\link{par}} to set graphical parameters before plotting.
    The \code{mfrow} splitting is handled per default. Afterwards,
    the \code{par}ameters are reverted to their original values.
    Use \code{par.list=NULL} to disable the internal \code{par} call.}

  \item{k}{the unit to plot, i.e., an element of \code{1:ncol(x)}.}
  \item{ylim}{the y limits of the plot(s). Ignored if
    \code{same.scale=FALSE}.}
  \item{axes}{a logical value indicating whether both axes should be drawn
    on the plot.}
  \item{xaxis.tickFreq,xaxis.labelFreq,xaxis.labelFormat}{see
    \code{\link{addFormattedXAxis}}.}
  \item{epochsAsDate}{Boolean indicating whether to treat the epochs as
    Date objects (or to transform them to dates such that the new x-axis
    formatting is applied).
    Default: Value of the \code{epochAsDate} slot of \code{x}.}
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
    elements to plot. See also the \code{col} argument.}
  \item{alarm.symbol}{a list with entries \code{pch}, \code{col},
    \code{cex} and \code{lwd} specifying the appearance of the outbreak
    symbol in the plot.}
  \item{outbreak.symbol}{a list with entries \code{pch}, \code{col},
    \code{cex} and \code{lwd} specifying the appearance of the outbreak
    symbol in the plot.}
  \item{legend.opts}{a list of arguments for \code{\link{legend}}.
    If \code{\link{missing}(legend.opts)} (i.e., not explicitly
    specified), the default legend will only be produced if \code{x}
    contains any information on outbreaks, alarms, or upperbounds.
    To disable the legend, use, e.g., \code{legend.opts=NULL}.
    Otherwise, the following arguments are default:
    \describe{
      \item{\code{x}}{\code{"top"}}
      \item{\code{legend}}{\code{c("Infected","Threshold","Outbreak","Alarm")}}
      \item{\code{lty,pch,col}}{the corresponding graphical settings}
    }
    Any further arguments to the \code{legend} function are
    just provided as additional elements of this list,
    e.g. \code{horiz=TRUE}.
  }
  \item{dx.upperbound}{horizontal change in the plotting of the
    upperbound line. Sometimes it can be convenient to offset this line
    a little for better visibility.}
  \item{lvl}{A vector of length \code{ncol(x)}, which is used to
    specify the hierarchy level for each time series in the sts object
    for alarm plots.}
  \item{cex}{A numerical value giving the amount by which plotting text
    and symbols should be magnified relative to the default. See
    \code{\link{par}} for details.}
  \item{cex.yaxis}{The magnification to be used for y-axis annotation
    relative to the current setting of \code{cex}.}
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
  Its arguments are (almost) similar to \code{\link{plot.survRes}}.
}

\value{
  \code{NULL} (invisibly).
  The functions are called for their side-effects.
}

\author{
  Michael H\enc{ö}{oe}hle and Sebastian Meyer
}

\seealso{
  There is an \code{\link[=autoplot.sts]{autoplot}}-method, which
  implements \CRANpkg{ggplot2}-based time-series plots of \code{"sts"}
  objects.
  
  The \code{\link{stsplot}} help page gives an overview of other
  types of plots for \code{"sts"} objects.
}

\examples{
data("ha.sts")
print(ha.sts)

plot(ha.sts, type=observed ~ time | unit)  # default multivariate type
plot(ha.sts, units=c("mitt", "pank"))      # selected units
plot(ha.sts, type=observed ~ time)         # aggregated over all districts

## Hook function example
hookFunc <- function() grid(NA,NULL,lwd=1)
plot(ha.sts, hookFunc=hookFunc)

## another multivariate time series example plotted "as.one"
data("measlesDE")
plot(measlesDE, units=1:2, as.one=TRUE, legend.opts=list(cex=0.8))
## more sophisticated plots are offered by package "xts"
if (require("xts"))
    plot(as.xts(measlesDE))

## Use ISO8601 date formatting (see ?strptime) and no legend
data("salmNewport")
plot(aggregate(salmNewport,by="unit"), xlab="Time (weeks)",
     xaxis.tickFreq=list("\%m"=atChange,"\%G"=atChange),
     xaxis.labelFreq=list("\%G"=atMedian),xaxis.labelFormat="\%G")

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
     xlab="Time (days)")
}

\keyword{hplot}
\keyword{ts}
