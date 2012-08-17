\encoding{latin1}
\name{twinstim_iafplot}
\alias{iafplot}

\title{
  Plot the spatial or temporal interaction function of a twimstim
}

\description{
  The function plots the fitted temporal or (isotropic) spatial
  interaction function of a \code{twinstim} object.
}

\usage{
iafplot(object, which = c("siaf", "tiaf"),
    types = 1:nrow(object$qmatrix),
    conf.type = if (length(pars) > 1) "bootstrap" else "parbounds",
    conf.level = 0.95, conf.B = 999,
    ngrid = 101, col.estimate = rainbow(length(types)), col.conf = col.estimate,
    alpha.B = 0.15, lwd = c(3,1), lty = c(1,2), xlim = c(0,eps), ylim = c(0,1),
    add = FALSE, xlab = NULL, ylab = NULL, ...)
}

\arguments{
  \item{object}{
    Object of class \code{"twinstim"} containing the fitted model.
  }
  \item{which}{
    Argument indicating which of the two interaction functions to plot.
    Possible values are \code{"siaf"} (default) for the spatial interaction
    function and \code{"tiaf"} for the temporal interaction function.
  }
  \item{types}{
    In case of a marked twinstim, for which event \code{types} 
    should the interaction function be plotted. 
  }
  \item{conf.type}{
    type of confidence interval to produce. If \code{"bootstrap"},
    \code{conf.B} parameter values (vectors) are sampled from the
    asymptotic (multivariate) normal distribution of the ML estimate(s) of the
    interaction function parameters; the interaction function is then
    evaluated on the grid of \code{ngrid} x-values (i.e. temporal or
    spatial distances from the host) for each parameter realization to
    obtain a \code{conf.level} confidence interval at each point of the x-grid.
    If \code{"parbounds"}, the \code{conf.level} Wald confidence
    intervals for the interaction function parameters are calculated and
    the interaction function is evaluated on the grid of \code{ngrid}
    x-values (distances from the host) for all combinations of the bounds
    of the parameters and the point-wise extremes of those functions are
    plotted. This type of confidence interval is actually only valid in
    the case of only one interaction function parameter, but could also
    be used as a rough indication if bootstrapping takes too long.
  }
  \item{conf.level}{
    the confidence level required. For \code{conf.type = "bootstrap"} it
    may also be specified as \code{NA}, in which case all \code{conf.B}
    bootstrapped functions will be plotted with transparency value given
    by \code{alpha.B}.
  }
  \item{conf.B}{
    number of samples for the \code{"bootstrap"} confidence interval.
  }
  \item{ngrid}{
    number of evaluation points of the interaction function in the
    interval \code{c(0,xlim[2])}. 
  }
  \item{col.estimate}{
    vector of colours to use for the function point estimates of the different \code{types}.
  }
  \item{col.conf}{
    vector of colours to use for the confidence intervals of the different \code{types}.
  }
  \item{alpha.B}{
    alpha transparency value (as relative opacity) used for the \code{conf.B}
    bootstrapped interaction functions in case \code{conf.level = NA}
  }
  \item{lwd, lty}{
    Numeric vectors of length two specifying the line width and type of point
    estimates (first element) and confidence limits (second element),
    respectively. 
  }
  \item{xlim, ylim}{
    Vectors of length two containing the x- and y-axis limit of the
    plot. The default y-axis ranges from 0 to 1.
    The default x-axis ranges from 0 to the length of the
    observation period (\code{which="tiaf"}) or the length of the
    diagonale of the bounding box of the observation region
    (\code{which="siaf"}), respectively.
  }
  \item{add}{
    Add to an existing plot?
  }
  \item{xlab, ylab}{
    labels for the axes with \code{NULL} providing sensible defaults.
  }
  \item{\dots}{
    Additional arguments passed to the default \code{plot} method.
  }
}

\value{
  A plot is created -- see e.g. Figure 3(b) in Meyer et al. (2012).
  \code{NULL} is returned (invisibly).
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{plot.twinstim}}, which calls this function.
}

\examples{
data(imdepifit)
iafplot(imdepifit, "tiaf", types=1)   # tiaf.constant(), not very exciting
iafplot(imdepifit, "siaf", xlim=c(0,200), types=1,   # same for types=2
     col.estimate=1, lwd=c(2,1))
# alternative: bootstrap-CI (actually only useful if iaf$npars > 1)
iafplot(imdepifit, "siaf", add=TRUE, types=1, col.estimate=0, ngrid=500,
     conf.type="bootstrap", col.conf=2)

\dontshow{
## FIXME: fit model with 2-parameter Lomax kernel as siaf (this takes a
##        while and does not converge properly...)
#data(imdepi)
#imdepifit_lomax <- update(imdepifit,
#    siaf = siaf.lomax(1),
#    optim.args = list(par=c(coef(imdepifit), 0), control=list(REPORT=1)),
#    nCub = 3, nCub.adaptive = FALSE, model = FALSE, cumCIF = FALSE
#)
}
}

\keyword{hplot}
\keyword{aplot}
