\encoding{latin1}
\name{twinstim_iafplot}
\alias{iafplot}

\title{
  Plot the Spatial or Temporal Interaction Function of a \code{twimstim}
}

\description{
  The function plots the fitted temporal or (isotropic) spatial
  interaction function of a \code{twinstim} object.
  The implementation is illustrated in Meyer et al. (2017, Section 3),
  see \code{vignette("twinstim")}.
}

\usage{
iafplot(object, which = c("siaf", "tiaf"), types = NULL,
        scaled = c("intercept", "standardized", "no"), truncated = FALSE,
        log = "", conf.type = if (length(pars) > 1) "MC" else "parbounds",
        conf.level = 0.95, conf.B = 999, xgrid = 101,
        col.estimate = rainbow(length(types)), col.conf = col.estimate,
        alpha.B = 0.15, lwd = c(3,1), lty = c(1,2),
        verticals = FALSE, do.points = FALSE,
        add = FALSE, xlim = NULL, ylim = NULL, xlab = NULL, ylab = NULL,
        legend = !add && (length(types) > 1), ...)
}

\arguments{
  \item{object}{
    object of class \code{"twinstim"} containing the fitted model.
  }
  \item{which}{
    argument indicating which of the two interaction functions to plot.
    Possible values are \code{"siaf"} (default) for the spatial interaction
    \eqn{f(x)} as a function of the distance \eqn{x}, and \code{"tiaf"}
    for the temporal interaction function \eqn{g(t)}.
  }
  \item{types}{
    integer vector indicating for which event \code{types} 
    the interaction function should be plotted in case of a marked
    \code{"twinstim"}. The default \code{types=NULL} checks if the interaction
    function is type-specific: if so, \code{types=1:nrow(object$qmatrix)}
    is used, otherwise \code{types=1}.
  }
  \item{scaled}{
    character string determining if/how the the interaction function
    should be scaled. Possible choices are:
    \describe{
      \item{"intercept":}{multiplication by the epidemic intercept.}
      \item{"standardized":}{division by the value at 0 distance such
        that the function starts at 1.}
      \item{"no":}{no scaling.}
    }
    The first one is the default and required for the comparison of
    estimated interaction functions from different models.
    For backward compatibility, \code{scaled} can also be a boolean,
    where \code{TRUE} refers to \code{"intercept"} scaling and
    \code{FALSE} to \code{"no"} scaling.
  }
  \item{truncated}{
    logical indicating if the plotted interaction function should
    take the maximum range of interaction (\code{eps.t}/\code{eps.s})
    into account, i.e., drop to zero at that point (if it is finite
    after all). If there is no common range of interaction, a
    \code{\link{rug}} indicating the various ranges will 
    be added to the plot if \code{truncated=TRUE}.
    If \code{truncated} is a scalar, this value is used as the point
    \code{eps} where the function drops to 0.
  }
  \item{log}{a character string passed to \code{\link{plot.default}}
    indicating which axes should be logarithmic.
    If \code{add=TRUE}, \code{log} is set according to
    \code{par("xlog")} and \code{par("ylog")}.}
  \item{conf.type}{
    type of confidence interval to produce.\cr
    If \code{conf.type="MC"} (or \code{"bootstrap"}), \code{conf.B}
    parameter vectors are sampled from the asymptotic
    (multivariate) normal distribution of the ML estimate of the
    interaction function parameters; the interaction function is then
    evaluated on the \code{xgrid} (i.e. temporal or spatial distances
    from the host) for each parameter realization to obtain a
    \code{conf.level} confidence interval at each point of the
    \code{xgrid} (or to plot the interaction functions of all
    Monte-Carlo samples if \code{conf.level=NA}).
    Note that the resulting plot is \code{\link{.Random.seed}}-dependent
    for the Monte-Carlo type of confidence interval.\cr
    If \code{conf.type="parbounds"}, the \code{conf.level} Wald confidence
    intervals for the interaction function parameters are calculated and
    the interaction function is evaluated on the \code{xgrid}
    (distances from the host) for all combinations of the bounds
    of the parameters and the point-wise extremes of those functions are
    plotted. This type of confidence interval is only valid in case of
    a single parameter, i.e. \code{scaled + nsiafpars == 1}, but could
    also be used as a rough indication if the Monte-Carlo approach takes
    too long. A warning is thrown if the \code{"parbounds"} type is used
    for multiple parameters.\cr 
    If \code{conf.type="none"} or \code{NA} or \code{NULL}, no
    confidence interval will be calculated.
  }
  \item{conf.level}{
    the confidence level required. For \code{conf.type = "MC"} it
    may also be specified as \code{NA}, in which case all \code{conf.B}
    sampled functions will be plotted with transparency value given
    by \code{alpha.B}.
  }
  \item{conf.B}{
    number of samples for the \code{"MC"} (Monte Carlo)
    confidence interval.
  }
  \item{xgrid}{
    either a numeric vector of x-values (distances from the host) where
    to evaluate \code{which}, or a scalar representing the desired number of
    evaluation points in the interval \code{c(0,xlim[2])}.\cr
    If the interaction function is a step function
    (\code{\link{siaf.step}} or \code{\link{tiaf.step}}), \code{xgrid}
    is ignored and internally set to \code{c(0, knots)}.
  }
  \item{col.estimate}{
    vector of colours to use for the function point estimates of the different \code{types}.
  }
  \item{col.conf}{
    vector of colours to use for the confidence intervals of the different \code{types}.
  }
  \item{alpha.B}{
    alpha transparency value (as relative opacity) used for the \code{conf.B}
    sampled interaction functions in case \code{conf.level = NA}
  }
  \item{lwd, lty}{
    numeric vectors of length two specifying the line width and type of point
    estimates (first element) and confidence limits (second element),
    respectively. 
  }
  \item{verticals,do.points}{graphical settings for step function
    kernels. These can be logical (as in \code{\link{plot.stepfun}}) or
    lists of graphical parameters.}
  \item{add}{
    add to an existing plot?
  }
  \item{xlim, ylim}{
    vectors of length two containing the x- and y-axis limit of the
    plot. The default y-axis range (\code{ylim=NULL}) is from 0 to the
    value of the (scaled) interaction function at \eqn{x = 0}.
    The default x-axis (\code{xlim=NULL}) starts at 0, and the upper
    limit is determined as follows (in decreasing order of precedence):
    \itemize{
      \item If \code{xgrid} is a vector of evaluation points, \code{xlim[2]} is
      set to \code{max(xgrid)}.
      \item \code{eps.t}/\code{eps.s} if it is unique and finite.
      \item If the interaction function is a step function with
      \code{maxRange<Inf}, i.e. it drops to 0 at \code{maxRange},
      \code{xlim[2]} is set to \code{maxRange}.
      \item Otherwise, it is set to the length of the observation period
      (\code{which="tiaf"}) or the diagonal length of the bounding box of
      the observation region (\code{which="siaf"}), respectively.
    }
  }
  \item{xlab, ylab}{
    labels for the axes with \code{NULL} providing sensible defaults.
  }
  \item{legend}{
    logical indicating if a legend for the \code{types} should be added.
    It can also be a list of arguments passed to \code{\link{legend}}
    to tweak the default settings.
  }
  \item{\dots}{
    additional arguments passed to the default \code{plot} method.
  }
}

\value{
  A plot is created -- see e.g. Figure 3(b) in Meyer et al. (2012).
  
  The function invisibly returns a matrix of the plotted values of the
  interaction function (evaluated on \code{xgrid}, by type). The first
  column of the matrix contains the distance \eqn{x}, and the remaining
  \code{length(types)} columns contain the (scaled) function values for
  each type.
  
  The pointwise confidence intervals of the interaction functions are
  returned in similar matrices as attributes: if
  \code{length(types)==1}, there is a single attribute \code{"CI"},
  whereas for multiple types, the attributes are named
  \code{paste0("CI.",typeNames)} (where the \code{typeNames} are
  retrieved from \code{object$qmatrix}).
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}

  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
}

\author{
  Sebastian Meyer
}

\seealso{
  \code{\link{plot.twinstim}}, which calls this function.
}

\examples{
data("imdepifit")

iafplot(imdepifit, "tiaf", scaled=FALSE)   # tiaf.constant(), not very exciting
iafplot(imdepifit, "siaf", scaled=FALSE)

# scaled version uses a Monte-Carlo-CI
set.seed(1)  # result depends on .Random.seed
iafplot(imdepifit, "siaf", scaled=TRUE, conf.type="MC", conf.B=199,
        col.conf=gray(0.4), conf.level=NA)  # show MC samples
}

\keyword{hplot}
\keyword{aplot}
