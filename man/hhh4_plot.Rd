\encoding{latin1}
\name{plot.hhh4}
\alias{plot.hhh4}
\alias{plotHHH4_fitted}
\alias{plotHHH4_fitted1}
\alias{plotHHH4_season}
\alias{getMaxEV_season}
\alias{plotHHH4_maxEV}
\alias{getMaxEV}
\alias{plotHHH4_maps}
\alias{plotHHH4_ri}
\alias{plotHHH4_neweights}

\title{Plots for Fitted \code{hhh4}-models}

\description{
  There are six \code{type}s of plots for fitted \code{\link{hhh4}} models:
  \itemize{
    \item Plot the \code{"fitted"} component means (of selected units)
    along time along with the observed counts.
    \item Plot the estimated \code{"season"}ality of the three components.
    \item Plot the time-course of the dominant eigenvalue \code{"maxEV"}.
    \item If the units of the corresponding multivariate
    \code{"\linkS4class{sts}"} object represent different regions,
    maps of the fitted mean components averaged over time (\code{"maps"}),
    or a map of estimated region-specific intercepts (\code{"ri"}) of a
    selected model component can be produced.
    \item Plot the (estimated) neighbourhood weights
    (\code{"neweights"}) as a function of neighbourhood order
    (shortest-path distance between regions), i.e., \code{w_ji ~ o_ji}.
  }
  Spatio-temporal \code{"hhh4"} models and these plots are illustrated in
  Meyer et al. (2017, Section 5), see \code{vignette("hhh4_spacetime")}.
}

\usage{
\method{plot}{hhh4}(x, type=c("fitted", "season", "maxEV", "maps", "ri", "neweights"), ...)

plotHHH4_fitted(x, units = 1, names = NULL,
                col = c("grey85", "blue", "orange"),
                pch = 19, pt.cex = 0.6, pt.col = 1,
                par.settings = list(),
                legend = TRUE, legend.args = list(),
                legend.observed = FALSE,
                decompose = NULL, total = FALSE, meanHHH = NULL, ...)

plotHHH4_fitted1(x, unit = 1, main = NULL,
                 col = c("grey85", "blue", "orange"),
                 pch = 19, pt.cex = 0.6, pt.col = 1, border = col, 
                 start = x$stsObj@start, end = NULL, xaxis = NULL,
                 xlim = NULL, ylim = NULL, xlab = "", ylab = "No. infected",
                 hide0s = FALSE, decompose = NULL, total = FALSE, meanHHH = NULL)

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

plotHHH4_maps(x, which = c("mean", "endemic", "epi.own", "epi.neighbours"),
              prop = FALSE, main = which, zmax = NULL, col.regions = NULL,
              labels = FALSE, sp.layout = NULL, ...,
              map = x$stsObj@map, meanHHH = NULL)

plotHHH4_ri(x, component, exp = FALSE,
            at = list(n = 10), col.regions = cm.colors(100),
            colorkey = TRUE, labels = FALSE, sp.layout = NULL,
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
    or \code{"maps"} of the fitted mean components averaged over time,
    or a map of estimated region-specific random
    intercepts (\code{"ri"}) of a specific model \code{component}.
    The latter two require \code{x$stsObj} to contain a map.}
  \item{\dots}{For \code{plotHHH4_season} and \code{plotHHH4_maxEV},
    one or more \code{\link{hhh4}}-fits, or a single list of these.
    Otherwise further arguments passed on to other functions.\cr
    For the \code{plot}-method these go to the specific plot
    \code{type} function.\cr
    \code{plotHHH4_fitted} passes them to \code{plotHHH4_fitted1}, which is called
    sequentially for every unit in \code{units}.\cr
    \code{plotHHH4_maps} and \code{plotHHH4_ri} pass additional arguments to
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
  \item{col,border}{length 3 vectors specifying the fill and border colors for the
    endemic, autoregressive, and spatio-temporal component polygons (in
    this order).}
  \item{pch,pt.cex,pt.col}{style specifications for the dots drawn to represent
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
  \item{decompose}{if \code{TRUE} or (a permutation of)
    \code{colnames(x$stsObj)}, the fitted mean will be decomposed
    into the contributions from each single unit and the endemic part
    instead of the default endemic + AR + neighbours decomposition.}
  \item{total}{logical indicating if the fitted components should be
    summed over all units to be compared with the total observed
    counts at each time point. If \code{total=TRUE}, the
    \code{units}/\code{unit} argument is ignored.}
  \item{start,end}{time range to plot specified by vectors of length two
    in the form \code{c(year,number)}, see \code{"\linkS4class{sts}"}.}
  \item{xaxis}{if this is a list (of arguments for
    \code{\link{addFormattedXAxis}}, the time axis is nicely labelled
    similar to \code{\link{stsplot_time}}.
    Note that in this case, the time indexes \code{1:nrow(x$stsObj)}
    will be used as x-values in the plot, which is different from the
    long-standing default (\code{xaxis = NULL}) with a real time scale.}
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
  \item{intercept}{logical indicating whether to include the global
    intercept. For \code{plotHHH4_season}, the default (\code{FALSE})
    means to plot seasonality as a multiplicative effect on the
    respective component. Multiplication by the intercept only
    makes sense if there are no further (non-centered)
    covariates/offsets in the component.}
  \item{exp}{logical indicating whether to \code{exp}-transform random
    effects to show multiplicative effects on the respective components.
    The default is \code{FALSE}.}
  \item{at}{a numeric vector of breaks for the color levels (see
    \code{\link[lattice]{levelplot}}), or a list specifying the number
    of breaks \code{n} (default: 10) and their \code{range} (default:
    range of the random effects, extended to be symmetric around 0,
    or around 1 if \code{exp=TRUE}). If \code{exp=TRUE}, the breaks are
    generated using \code{scales::\link[scales]{log_breaks}}.}
  \item{matplot.args}{list of line style specifications passed to
    \code{\link{matplot}}, e.g., \code{lty}, \code{lwd}, \code{col}.}
  \item{refline.args}{list of line style specifications (e.g.,
    \code{lty} or \code{col}) passed to \code{\link{abline}} when
    drawing the reference line (\code{h=1}) in plots of seasonal effects
    (if \code{intercept=FALSE}) and of the dominant eigenvalue.
    The reference line is omitted if \code{refline.args} is not a list.}
  \item{which}{a character vector specifying the components of the mean
    for which to produce maps. By default, the overall mean and all
    three components are shown.}
  \item{prop}{a logical indicating whether the component maps should
    display proportions of the total mean instead of absolute numbers.}
  \item{zmax}{a numeric vector of length \code{length(which)} (recycled
    as necessary) specifying upper limits for the color keys of the
    maps, using a lower limit of 0.
    A missing element (\code{NA}) means to use a map-specific color key
    only covering the range of the values in that map
    (can be useful for \code{prop = TRUE}).
    The default \code{zmax = NULL} means to use the same scale for the component maps
    and a separate scale for the map showing the overall mean.}
  \item{col.regions}{a vector of colors used to encode the fitted
    component means (see \code{\link[lattice]{levelplot}}).
    For \code{plotHHH4_maps}, the length of this color vector also
    determines the number of levels, using 10 heat colors by default.}
  \item{colorkey}{a Boolean indicating whether to draw the color key.
    Alternatively, a list specifying how to draw it, see
    \code{\link[lattice]{levelplot}}.}
  \item{map}{an object inheriting from \code{"\linkS4class{SpatialPolygons}"}
    with \code{row.names} covering \code{colnames(x)}.}
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
  \code{plotHHH4_maps} returns a \code{\link[lattice]{trellis.object}} if
  \code{length(which) == 1} (a single \code{\link{spplot}}), and
  otherwise uses \code{\link[gridExtra]{grid.arrange}} from the
  \pkg{gridExtra} package to arrange all \code{length(which)}
  \code{\link{spplot}}s on a single page.
  \code{plotHHH4_ri} returns the generated \code{\link{spplot}}, i.e.,
  a \code{\link[lattice]{trellis.object}}.\cr
  \code{plotHHH4_neweights} eventually calls \code{plotter} and
  thus returns whatever is returned by that function.
}

\author{
  Sebastian Meyer
}

\references{
  Held, L. and Paul, M. (2012): Modeling seasonality in space-time
  infectious disease surveillance data.
  \emph{Biometrical Journal}, \bold{54}, 824-843.
  \doi{10.1002/bimj.201200037}

  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
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

## sum fitted components over all units
plot(measlesFit, total=TRUE)

## 'xaxis' option for a nicely formatted time axis
## default tick locations and labels:
plot(measlesFit, total=TRUE, xaxis=list(epochsAsDate=TRUE, line=1))
## an alternative with monthly ticks:
oopts <- surveillance.options(stsTickFactors = c("\%m"=0.75, "\%Y" = 1.5))
plot(measlesFit, total=TRUE, xaxis=list(epochsAsDate=TRUE,
    xaxis.tickFreq=list("\%m"=atChange, "\%Y"=atChange),
    xaxis.labelFreq=list("\%Y"=atMedian), xaxis.labelFormat="\%Y"))
surveillance.options(oopts)

## plot the multiplicative effect of seasonality
plot(measlesFit, type="season")

## dominant eigenvalue of the Lambda matrix (cf. Held and Paul, 2012)
getMaxEV(measlesFit)  # here simply constant and equal to exp(ar.1)
plot(measlesFit, type="maxEV")  # not very exciting

## fitted mean components/proportions by district, averaged over time
if (requireNamespace("gridExtra")) {
    plot(measlesFit, type="maps", labels=list(cex=0.6),
         which=c("mean", "endemic", "epi.own"),
         main=c("total", "endemic", "autoregressive"))
    plot(measlesFit, type="maps", labels=list(cex=0.6),
         which=c("endemic", "epi.own"), prop=TRUE, zmax=NA,
         main=c("endemic proportion", "autoregressive proportion"))
}

## estimated random intercepts of the endemic component
fixef(measlesFit)["end.ri(iid)"]     # global intercept (log-scale)
ranef(measlesFit, tomatrix = TRUE)   # zero-mean deviations
ranef(measlesFit, intercept = TRUE)  # sum of the above
exp(ranef(measlesFit))               # multiplicative effects
plot(measlesFit, type="ri", component="end",
     main="deviations around the endemic intercept (log-scale)")
plot(measlesFit, type="ri", component="end", exp=TRUE,
     main="multiplicative effects", labels=list(font=3, labels="GEN"))

## neighbourhood weights as a function of neighbourhood order
plot(measlesFit, type="neweights")  # boring, model has no "ne" component

## fitted values for the 6 regions with most cases and some customization
bigunits <- tail(names(sort(colSums(observed(measlesWeserEms)))), 6)
plot(measlesFit, units=bigunits,
     names=measlesWeserEms@map@data[bigunits,"GEN"],
     legend=5, legend.args=list(x="top"), xlab="Time (weekly)",
     hide0s=TRUE, ylim=c(0,max(observed(measlesWeserEms)[,bigunits])),
     start=c(2002,1), end=c(2002,26), par.settings=list(xaxs="i"))

## plot completely decomposed mean structure (useless without 'ne' component)
plot(measlesFit, units=bigunits, col=rainbow(measlesFit$nUnit), decompose=TRUE)
}

\keyword{hplot}
