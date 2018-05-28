\name{plot.survRes}
\alias{plot.survRes}
\alias{plot.survRes.one}
\encoding{latin1}

\title{Plot a survRes object}
\description{
  Plotting a (multivariate) \code{survRes} object. The internal
  function \code{plot.survRes.one} is used as a helper function
  to plot a univariate time series.
}
\usage{
  \method{plot}{survRes}(x, method=x$control$name, disease=x$control$data, 
  xaxis.years=TRUE,startyear = 2001, firstweek = 1, same.scale=TRUE, ...)

  plot.survRes.one(x, method=x$control$name, disease=x$control$data,
   domany=FALSE,ylim=NULL,xaxis.years=TRUE,startyear = 2001, firstweek = 1, 
   xlab="time", ylab="No. infected", main=NULL, type="hhs",
   lty=c(1,1,2),col=c(1,1,4), 
   outbreak.symbol = list(pch=3,col=3),alarm.symbol=list(pch=24,col=2),
   legend.opts=list(x="top",
   legend=c("Infected", "Upperbound", "Alarm", "Outbreak"),
   lty=NULL,col=NULL,pch=NULL), ...)
}
\arguments{
        \item{x}{object of class \code{survRes}}
        \item{method}{surveillance method to be used in title}
        \item{disease}{name of disease in title}
        \item{xaxis.years}{Boolean indicating whether to show a year based x-axis for weekly data} 
        \item{domany}{Boolean telling the function whether it is called for a multivariate (\code{TRUE}) or univariate (\code{FALSE}) \code{survRes} object. In case of \code{TRUE} no titles are drawn.}
        \item{ylim}{range of y axis}
        \item{startyear}{year to begin the axis labeling (the year where
                the oldest data come from)}
        \item{firstweek}{number of the first week of January in the first year (just for axis labeling reasons)}
        \item{xlab}{label of the x-axis}
        \item{ylab}{label of the y-axis}
        \item{main}{the title of the graphics is generated from the \code{method} and \code{disease} arguments if not specified otherwise}
        \item{same.scale}{plot all time series with the same \code{ylim}? Defaults to \code{true}}.
        \item{type}{line type of the observed counts (first two elements) and the upper bound (third element)}
        \item{lty}{vector of size 3 specifying the line type of the observed counts (left, right) and the upperbound line}
        \item{col}{vector with three elements: color of left bar and color of top bar, color of right bar, col of the upperbound line.}
        \item{outbreak.symbol}{list with entries \code{pch} and \code{col} specifying the plot symbol}
        \item{alarm.symbol}{list with entries \code{pch} and \code{col} specifying the plot symbol}
        \item{legend.opts}{a list containing the entries to be sent to the \code{\link{legend}} function. If no legend is requested use \code{legend.opts=NULL}. Otherwise, the following arguments are default
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
        \item{...}{arguments passed to \code{plot.survRes.one}.
            From there, further arguments are passed to \code{\link{matplot}}.}
}
\value{
  none. A plot showing the number of infected, the threshold for recognizing an
  outbreak, the alarm status and the outbreak status is generated.
}
\author{M. \enc{Höhle}{Hoehle}}
\examples{
data(ha)
ctrl <- list(range = 209:290, b = 2, w = 6, alpha = 0.005)
plot(algo.bayes(aggregate(ha), control = ctrl))
}

\keyword{hplot}
