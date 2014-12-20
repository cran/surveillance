\name{plot.hhh4}
\alias{plot.hhh4}
\alias{plotHHH4_fitted}
\alias{plotHHH4_fitted1}
\alias{plotHHH4_season}
\alias{getMaxEV_season}
\alias{plotHHH4_maxEV}
\alias{getMaxEV}
\alias{plotHHH4_ri}
\alias{plotHHH4_neweights}

\title{Plots for Fitted \code{hhh4}-models}

\description{
  There are five \code{type}s of plots for fitted \code{\link{hhh4}} models:
  \itemize{
    \item Plot the \code{"fitted"} component means (of selected units)
    along time along with the observed counts.
    \item Plot the estimated \code{"season"}ality of the three components.
    \item Plot the time-course of the dominant eigenvalue \code{"maxEV"}.
    \item If the units of the corresponding multivariate
    \code{"\linkS4class{sts}"} object represent different regions, a map
    of estimated region-specific intercepts (\code{"ri"}) of a selected
    model component can be produced.
    \item Plot the (estimated) neighbourhood weights
    (\code{"neweights"}) as a function of neighbourhood order
    (shortest-path distance between regions), i.e., \code{w_ji ~ o_ji}.
  }
}

\usage{
\method{plot}{hhh4}(x, type=c("fitted", "season", "maxEV", "ri", "neweights"), ...)

plotHHH4_fitted(x, units = 1, names = NULL,
                col = c("orange", "blue", "grey85", "black"),
                pch = 19, pt.cex = 0.6,
                par.settings = list(),
                legend = TRUE, legend.args = list(),
                legend.observed = FALSE, ...)

plotHHH4_fitted1(x, unit = 1, main = NULL,
                 col = c("grey30", "grey60", "grey85", "grey0"),
                 pch = 19, pt.cex = 0.6, border = col, 
                 start = x$stsObj@start, end = NULL,
                 xlim = NULL, ylim = NULL, xlab = "", ylab = "No. infected",
                 hide0s = FALSE, meanHHH = NULL)

plotHHH4_season(..., components = NULL, intercept = FALSE,
                xlim = NULL, ylim = NULL,
                xlab = NULL, ylab = "", main = NULL,
                par.settings = list(), matplot.args = list(),
                legend = NULL, legend.args = list(),
                refline.args = list(), unit = 1)
getMaxEV_season(x)

plotHHH4_maxEV(...,
               matplot.args = list(), refline.args = list(),
               legend.args = list())
getMaxEV(x)

plotHHH4_ri(x, component, labels = FALSE, sp.layout = NULL,
            gpar.missing = list(col = "darkgrey", lty = 2, lwd = 2),
            ...)

plotHHH4_neweights(x, plotter = boxplot, ...,
                   exclude = 0, maxlag = Inf)
}

\arguments{
  \item{x}{a fitted \code{\link{hhh4}} object.}
  \item{type}{type of plot: either \code{"fitted"} component means of
    selected \code{units} along time along with the observed counts, or
    \code{"season"}ality plots of the model components and the epidemic
    dominant eigenvalue (which may also be plotted along overall time by
    \code{type="maxEV"}, especially if the model contains time-varying
    neighbourhood weights or unit-specific epidemic effects),
    or a map of estimated region-specific random
    intercepts (\code{"ri"}) of a specific model \code{component}.
    The latter requires the model \code{x} to contain the selected
    random intercepts, and \code{x$stsObj} to contain a map.}
  \item{\dots}{For \code{plotHHH4_season} and \code{plotHHH4_maxEV},
    one or more \code{\link{hhh4}}-fits, or a single list of these.
    Otherwise further arguments passed on to other functions.\cr
    For the \code{plot}-method these go to the specific plot
    \code{type} function.\cr
    \code{plotHHH4_fitted} passes them to \code{plotHHH4_fitted1}, which is called
    sequentially for every unit in \code{units}.\cr
    \code{plotHHH4_ri} passes additional arguments to
    \code{\link{spplot}}, and \code{plotHHH4_neweights} to the
    \code{plotter}.}
  \item{units,unit}{integer or character vector specifying a single
    \code{unit} or possibly multiple \code{units} to plot.
    It indexes \code{colnames(x$stsObj)}.\cr
    In \code{plotHHH4_fitted}, \code{units=NULL} plots all units.\cr
    In the seasonality plot,
    selection of a unit is only relevant if the model contains
    unit-specific intercepts or seasonality terms.}
  \item{names,main}{main title(s) for the selected
    \code{unit}(\code{s}) / \code{components}. If \code{NULL} (default),
    \code{plotHHH4_fitted1} will use the appropriate element of
    \code{colnames(x$stsObj)}, whereas \code{plotHHH4_season} uses
    default titles.}
  \item{col,border}{length 4 vector specifying the colors for the
    spatio-temporal, autoregressive, endemic, and observed elements in
    this order. For \code{plotHHH4_fitted}, a length 3 vector also
    works, in which case the default \code{"black"} will be used for the
    dots of observed counts.}
  \item{pch,pt.cex}{style specifications for the dots drawn to represent
    the observed counts. \code{pch=NA} can be used to disable these dots.}
  \item{par.settings}{list of graphical parameters for
    \code{\link{par}}. Sensible defaults for \code{mfrow}, \code{mar} and
    \code{las} will be applied unless overridden or
    \code{!is.list(par.settings)}.}
  \item{legend}{Integer vector specifying in which of the
    \code{length(units)} frames the legend should be drawn. If a logical
    vector is supplied, \code{which(legend)} determines the frame
    selection, i.e., the default is to drawn the legend in the first
    (upper left) frame only, and \code{legend=FALSE} results in no
    legend being drawn.}
  \item{legend.args}{list of arguments for \code{\link{legend}}, e.g.,
    to modify the default positioning
    \code{list(x="topright", inset=0.02)}.}
  \item{legend.observed}{logical indicating if the legend should contain
    a line for the dots corresponding to observed counts.}
  \item{start,end}{time range to plot specified by vectors of length two
    in the form \code{c(year,number)}, see \code{"\linkS4class{sts}"}.}
  \item{xlim}{numeric vector of length 2 specifying the x-axis range.
    The default (\code{NULL}) is to plot the complete time range.}
  \item{ylim}{y-axis range.
    For \code{type="fitted"}, this defaults to 
    \code{c(0,max(observed(x$stsObj)[,unit]))}.
    For \code{type="season"}, \code{ylim} must be a list of length
    \code{length(components)} specifying the range for every component
    plot, or a named list to customize only a subset of these. If only
    one \code{ylim} is specified, it will be recycled for all
    \code{components} plots.}
  \item{xlab,ylab}{axis labels. For \code{plotHHH4_season}, \code{ylab}
    specifies the y-axis labels for all \code{components} in a
    list (similar to \code{ylim}). If \code{NULL} or incomplete, 
    default mathematical expressions are used.
    If a single name is supplied such as the default \code{ylab=""} (to
    omit y-axis labels), it is used for all \code{components}.}
  \item{hide0s}{logical indicating if dots for zero observed counts
    should be omitted. Especially useful if there are too many.}
  \item{meanHHH}{(internal) use different component means than those
    estimated and available from \code{x}.}
  \item{components}{character vector of component names, i.e., a subset
    of \code{c("ar", "ne", "end")}, for which to plot the estimated
    seasonality. If \code{NULL} (the default), only components which
    appear in any of the models in \code{\dots} are plotted.\cr
    A seasonality plot of the epidemic dominant eigenvalue
    is also available by including \code{"maxEV"} in \code{components},
    but it only supports models without epidemic covariates/offsets.}
  \item{intercept}{logical indicating whether to plot seasonality as a
    multiplicative effect on the respective component (the default
    \code{intercept=FALSE}), or additionally multiplied by the
    corresponding intercept (\code{intercept=TRUE}). The latter only
    makes sense if there are no further (non-centered)
    covariates/offsets in the component.}
  \item{matplot.args}{list of line style specifications passed to
    \code{\link{matplot}}, e.g., \code{lty}, \code{lwd}, \code{col}.}
  \item{refline.args}{list of line style specifications passed to
    \code{\link{abline}} to draw the reference line in the plot of
    the dominant eigenvalue, e.g., \code{lty} and \code{col}.}
  \item{component}{component for which to plot the estimated
    region-specific random intercepts. Must partially match one of
    \code{colnames(ranef(x, tomatrix=TRUE))}.}
  \item{labels}{determines if and how regions are labeled, see
    \code{\link{layout.labels}}.}
  \item{sp.layout}{optional list of additional layout items, see
    \code{\link{spplot}}.}
  \item{gpar.missing}{list of graphical parameters for
    \code{\link{sp.polygons}}, applied to regions
    with missing random intercepts, i.e., not included in the model.
    Such extra regions won't be plotted if
    \code{!is.list(gpar.missing)}.}
  \item{plotter}{the (name of a) function used to produce the plot of
    weights (a numeric vector) as a function of neighbourhood order (a
    factor variable). It is called as
    \code{plotter(Weight ~ Distance, ...)} and defaults to
    \code{\link{boxplot}}. A useful alternative is, e.g.,
    \code{\link{stripplot}} from package \pkg{lattice}.}
  \item{exclude}{vector of neighbourhood orders to be excluded from
    plotting (passed to \code{\link{factor}}).
    By default, the neighbourhood weight for order 0 is not
    shown, which is usually zero anyway.}
  \item{maxlag}{maximum order of neighbourhood to be assumed when
    computing the \code{\link{nbOrder}} matrix. This additional step is
    necessary iff \code{neighbourhood(x$stsObj)} only specifies a binary
    adjacency matrix.}
}

\value{
  \code{plotHHH4_fitted1} invisibly returns a matrix of the fitted
  component means for the selected \code{unit}, and \code{plotHHH4_fitted}
  returns these in a list for all \code{units}.\cr
  \code{plotHHH4_season} invisibly returns the plotted y-values, i.e. the
  multiplicative seasonality effect within each of \code{components}.
  Note that this will include the intercept, i.e. the point estimate of
  \eqn{exp(intercept + seasonality)} is plotted and returned.\cr
  \code{getMaxEV_season} returns a list with elements
  \code{"maxEV.season"} (as plotted by
  \code{plotHHH4_season(..., components="maxEV")},
  \code{"maxEV.const"} and \code{"Lambda.const"} (the Lambda matrix and
  its dominant eigenvalue if time effects are ignored).\cr
  \code{plotHHH4_maxEV} (invisibly) and \code{getMaxEV} return the
  dominant eigenvalue of the \eqn{\Lambda_t} matrix for all time points
  \eqn{t} of \code{x$stsObj}.\cr
  \code{plotHHH4_ri} returns the generated \code{\link{spplot}}, which
  is a \pkg{lattice} plot of class
  \code{"\link[=trellis.object]{trellis}"}.\cr
  \code{plotHHH4_neweights} eventually calls \code{plotter} and
  thus returns whatever is returned by that function.
}

\author{
  Sebastian Meyer
}

\references{
  Held, L. and Paul, M. (2012): Modeling seasonality in space-time
  infectious disease surveillance data.
  \emph{Biometrical Journal}, \bold{54}, 824-843.\cr
  DOI-Link: \url{http://dx.doi.org/10.1002/bimj.201200037}
}

\seealso{
  other methods for \code{hhh4} fits, e.g., \code{\link{summary.hhh4}}.
}

\examples{
data("measlesWeserEms")

## fit a simple hhh4 model
measlesModel <- list(
    ar = list(f = ~ 1),
    end = list(f = addSeason2formula(~0 + ri(type="iid"), S=1, period=52),
               offset = population(measlesWeserEms)),
    family = "NegBin1"
    )
measlesFit <- hhh4(measlesWeserEms, measlesModel)

## fitted values for a single unit
plot(measlesFit, units=2)

## plot the multiplicative effect of seasonality
plot(measlesFit, type="season")

## dominant eigenvalue of the Lambda matrix (cf. Held and Paul, 2012)
getMaxEV(measlesFit)  # here simply constant and equal to exp(ar.1)
plot(measlesFit, type="maxEV")  # not very exciting

## random intercepts of the endemic component
plot(measlesFit, type="ri", component="end", labels=list(font=3, labels="GEN"))

## neighbourhood weights as a function of neighbourhood order
plot(measlesFit, type="neweights")  # boring, model has no "ne" component

## fitted values for the 6 regions with most cases and some customization
bigunits <- tail(names(sort(colSums(observed(measlesWeserEms)))), 6)
plot(measlesFit, units=bigunits,
     names=measlesWeserEms@map@data[bigunits,"GEN"],
     legend=5, legend.args=list(x="top"), xlab="Time (weekly)",
     hide0s=TRUE, ylim=c(0,max(observed(measlesWeserEms)[,bigunits])),
     start=c(2002,1), end=c(2002,26), par.settings=list(xaxs="i"))
}

\keyword{hplot}
