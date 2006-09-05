\name{plot.survRes}
\alias{plot.survRes}
\alias{plot.survRes.one}
\title{Plot a survRes object}
\description{Plotting of a (multivariate) \code{survRes} object. The function \code{plot.survRes.one} is used as a helper function to plot a univariate time series.}
}
\usage{
    plot.survRes(x, method=x$control$name, disease=x$control$data, xaxis.years=TRUE,startyear = 2001, firstweek = 1, legend=TRUE,same.scale=TRUE,...)
    plot.survRes.one(x, method=x$control$name, disease=x$control$data, domany=FALSE,ylim=NULL,xaxis.years=TRUE,startyear = 2001, firstweek = 1, legend=TRUE, xlegpos =1/4, ylegpos = 1,...)

}
\arguments{
        \item{x}{Object of class \code{survRes}}
        \item{method}{Surveillance method to be used in title}
        \item{disease}{Name of disease in title}
        \item{xaxis.years}{Boolean indicating whether to show a year based x-axis for weekly data} 
          \item{domany}{Bollean telling the function is called for a multivariate (\code{true}) or univariate (\code{false}) \code{survRes} object.}
        \item{startyear}{Year to begin the axis labeling (the year where
                the oldest data come from)}
        \item{firstweek}{Number of the first week of January in the first year (just for axis labeling reasons)}
        \item{legend}{Boolean indicating whether to add a legend}
        \item{same.scale}{Plot all time series with the same \code{ylim}? Defaults to \code{true}}.
        \item{ylim}{range of y axis}
        \item{xlegpos}{x position of the legend as fraction of \code{max(xlim)}}
        \item{ylegpos}{y position of legend as fraction of \code{max(ylim)}}
        \item{...}{further arguments for the function \code{matplot}}
}
\value{
        \item{a plot}{showing the number of infected, the threshold for recognizing an
                outbreak, the computed alarm status and the defined alarm status.
        }
}
\details{
  The \code{plot.survRes.one} is intended for internal use. At the moment none of the surveillance methods support multivariate \code{survRes} objects. New versions of the packages currently under development will handle this.
}
\author{M. Höhle}
\examples{
data(ha)
ctrl <- list(range = 209:290, b = 2, w = 6, alpha = 0.005)
plot(algo.bayes(aggregate(ha), control = ctrl))
}

\keyword{}














% magic.dim
