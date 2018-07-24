\name{fanplot}
\alias{fanplot}
\title{Fan Plot of Forecast Distributions}

\description{
The \code{fanplot()} function in \pkg{surveillance} wraps functionality of
the dedicated \CRANpkg{fanplot} package, employing a different default style
and optionally adding point predictions and observed values.
}

\usage{
fanplot(quantiles, probs, means = NULL, observed = NULL, start = 1,
        fan.args = list(), means.args = list(), observed.args = list(),
        key.args = NULL, xlim = NULL, ylim = NULL, log = "",
        xlab = "Time", ylab = "No. infected", add = FALSE, ...)
}

\arguments{
  \item{quantiles}{
    a time x \code{probs} matrix of forecast quantiles at each time point.
  }
  \item{probs}{
    numeric vector of probabilities with values between 0 and 1.
  }
  \item{means}{
    (optional) numeric vector of point forecasts.
  }
  \item{observed}{
    (optional) numeric vector of observed values.
  }
  \item{start}{
    time index (x-coordinate) of the first prediction.
  }
  \item{fan.args}{
    a list of graphical parameters for the \code{\link[fanplot]{fan}},
    e.g., to employ a different \code{\link{colorRampPalette}} as
    \code{fan.col}, or to enable contour lines via \code{ln}.
  }
  \item{means.args}{
    a list of graphical parameters for \code{\link{lines}}
    to modify the plotting style of the \code{means}.
    The default is a white line within the fan.
  }
  \item{observed.args}{
    a list of graphical parameters for \code{\link{lines}}
    to modify the plotting style of the \code{observed} values.
  }
  \item{key.args}{
    if a list, a color key (in \code{\link[fanplot]{fan}()}'s
    \code{"boxfan"}-style) is added to the fan chart. The list may
    include positioning parameters \code{start} (the x-position) and
    \code{ylim} (the y-range of the color key), \code{space} to modify
    the width of the boxfan, and \code{rlab} to modify the labels.
    An alternative way of labeling the quantiles is
    via the argument \code{ln} in \code{fan.args}.
  }
  \item{xlim,ylim}{
    axis ranges.
  }
  \item{log}{
    a character string specifying which axes are to be logarithmic,
    e.g., \code{log="y"} (see \code{\link{plot.default}}).
  }
  \item{xlab,ylab}{
    axis labels.
  }
  \item{add}{
    logical indicating if the fan plot should be added to an existing plot.
  }
  \item{\dots}{
    further arguments are passed to \code{\link{plot.default}}.
    For instance, \code{panel.first} could be used to initialize the
    plot with \code{\link{grid}(nx=NA, ny=NULL)} lines.
  }
}

\value{
  \code{NULL} (invisibly), with the side effect of drawing a fan chart.
}

\author{
  Sebastian Meyer
}

\seealso{
  the underlying \code{\link[fanplot]{fan}} function in package
  \CRANpkg{fanplot}.
  The function is used in \code{\link{plot.oneStepAhead}} and
  \code{\link{plot.hhh4sims}}.
}

\examples{
## artificial data example to illustrate the graphical options
if (requireNamespace("fanplot")) {
    means <- c(18, 19, 20, 25, 26, 35, 34, 25, 19)
    y <- rlnorm(length(means), log(means), 0.5)
    quantiles <- sapply(1:99/100, qlnorm, log(means), seq(.5,.8,length.out=length(means)))
    
    ## default style with point predictions, color key and log-scale
    fanplot(quantiles = quantiles, probs = 1:99/100, means = means,
            observed = y, key.args = list(start = 1, space = .3), log = "y")
    
    ## with contour lines instead of a key, and different colors
    pal <- colorRampPalette(c("darkgreen", "gray93"))
    fanplot(quantiles = quantiles, probs = 1:99/100, observed = y,
            fan.args = list(fan.col = pal, ln = c(5,10,25,50,75,90,95)/100),
            observed.args = list(type = "b", pch = 19))
}
}

\keyword{hplot}
\keyword{distribution}
