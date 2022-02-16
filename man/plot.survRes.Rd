\name{plot.survRes}
\alias{plot.survRes}
\title{Plot a \code{survRes} object}
\description{
  Plotting a (multivariate) \code{survRes} object (soft-deprecated),
  including the number of infected, the threshold for recognizing an
  outbreak, the alarm status and the outbreak status.

  As of \pkg{surveillance} 1.20.0, legacy \code{survRes} objects are
  plotted via internal \code{\link{disProg2sts}} conversion
  and \code{\link{stsplot_time}}.
}
\usage{
  \method{plot}{survRes}(x, method = x$control$name, disease = x$control$data,
     xaxis.years = TRUE, startyear = 2001, firstweek = 1,
     same.scale = TRUE, ...,
     main = paste0("Analysis of ", disease, " using ", method))
}
\arguments{
  \item{x}{object of class \code{survRes}}
  \item{method}{surveillance method to be used in title}
  \item{disease}{name of disease in title}
  \item{xaxis.years}{Boolean indicating whether to show a year-based x-axis}
  \item{startyear,firstweek}{(legacy arguments, ignored with a warning)}
  \item{same.scale}{if \code{TRUE} (default), all plots use the same \code{ylim}}
  \item{\dots}{further arguments passed to \code{\link{stsplot_time}}}
  \item{main}{the plot title is generated from the
    \code{method} and \code{disease} arguments if not specified otherwise}
}
\examples{
data(ha)
ctrl <- list(range = 209:290, b = 2, w = 6, alpha = 0.005)
plot(algo.bayes(aggregate(ha), control = ctrl))
}
\keyword{hplot}
\keyword{internal}
