\name{stsplot_space}
\alias{stsplot_space}
\title{
  Map of Disease Counts/Incidence accumulated over a Given Period
}

\description{
This is the \code{plot} variant of \code{type=observed~unit} for
\code{"\linkS4class{sts}"} objects, i.e.,
\code{plot(stsObj, type=observed~unit, ...)} calls the function
documented below. It produces an \code{\link[sp]{spplot}}
where regions are color-coded according to disease incidence
(either absolute counts or relative to population) over a given
time period.
}

\usage{
stsplot_space(x, tps = NULL, map = x@map, population = NULL,
              main = NULL, labels = FALSE, ...,
              at = 10, col.regions = NULL,
              colorkey = list(space = "bottom", labels = list(at=at)),
              total.args = NULL,
              gpar.missing = list(col = "darkgrey", lty = 2, lwd = 2),
              sp.layout = NULL, %aspect = mapasp(map),
              xlim = bbox(map)[1, ], ylim = bbox(map)[2, ])
}

\arguments{
  \item{x}{
    an object of class \code{"\linkS4class{sts}"} or a matrix of counts,
    i.e., \code{observed(stsObj)}, where especially
    \code{colnames(x)} have to be contained in \code{row.names(map)}.
    If a matrix, the \code{map} object has to be provided explicitly.
    The possibility of specifying a matrix is, e.g., useful to plot mean
    counts of simulations from \code{\link{simulate.hhh4}}.
  }
  \item{tps}{
    a numeric vector of one or more time points.
    The unit-specific \emph{sum} over all time points \code{tps} is
    plotted. The default \code{tps=NULL} means accumulation over the whole
    time period \code{1:nrow(x)}.
  }
  \item{map}{
    an object inheriting from \code{"\linkSPclass{SpatialPolygons}"}
    representing the \code{ncol(x)} regions. By default
    the \code{map} slot of \code{x} is queried (which might be
    empty and is not applicable if \code{x} is a matrix of counts).
  }
  \item{population}{
    if \code{NULL} (default), the map shows the region-specific numbers
    of cases accumulated over \code{tps}. For a disease incidence map,
    \code{population} can be specified in three ways:
    \itemize{
      \item a numeric vector of population numbers in the
      \code{ncol(x)} regions, used to divide the disease counts.
      \item a matrix of population counts of dimension \code{dim(x)}
      (such as \code{population(x)} in an \code{"sts"} object).
      This will produce the cumulative incidence over \code{tps}
      relative to the population at the first time point, i.e., only
      \code{population[tps[1],]} is used.
      \item [if \code{is(x, "sts")}]
      a scalar specifying how \code{population(x)} should be scaled for
      use as the population matrix, i.e.,
      \code{population(x)/population} is used. For instance, if
      \code{population(x)} contains raw population numbers,
      \code{population=1000} would produce the incidence per 1000
      inhabitants.
    }
  }
  \item{main}{
    a main title for the plot. If \code{NULL} and \code{x} is of
    class \code{"sts"}, the time range of \code{tps} is put as
    the main title.
  }
  \item{labels}{
    determines if and how the regions of the \code{map} are labeled,
    see \code{\link{layout.labels}}.
  }
  \item{\dots}{
    further arguments for \code{\link[sp]{spplot}},
    for example \code{col = "white"} for white polygon lines.
  }
  \item{at}{
    either a number of levels (default: 10) for the categorization
    (color-coding) of counts/incidence,
    or a numeric vector of specific break points,
    or a named list of a number of levels (\code{"n"}), a transformer
    (\code{"trafo"}) of class \code{"\link[scales]{trans}"} defined by
    package \pkg{scales}, and optional further arguments for
    \code{\link{pretty}}. The default breaks are equally spaced
    on the square-root scale (equivalent to
    \code{\link[scales]{sqrt_trans}}).
    Note that intervals given by \code{at} are closed on the left
    and open to the right; if manually specified break points do not
    cover the data range, further breaks are automatically added at 0
    and the maximum (rounded up to 1 significant digit), respectively.
  }
  \item{col.regions}{
    a vector of fill colors, sufficiently long to serve all levels
    (determined by \code{at}).
    \dQuote{Heat} colors are used by default (\code{NULL}).
  }
  \item{colorkey}{
    a list describing the color key, see
    \code{\link[lattice]{levelplot}}. The default list elements will be
    updated by the provided list using \code{\link{modifyList}}.
  }
  \item{total.args}{
    an optional list of arguments for \code{\link[grid]{grid.text}} to
    have the overall number/incidence of cases printed at an edge of the
    map. The default settings are
    \code{list(label="Overall: ", x=1, y=0)}, and
    \code{total.args=list()} will use all of them.
  }
  \item{gpar.missing}{list of graphical parameters for
    \code{\link[sp]{sp.polygons}} applied to the regions of \code{map},
    which are not part of \code{x}. Such extra regions won't be plotted
    if \code{!is.list(gpar.missing)}.}
  \item{sp.layout}{
    optional list of additional layout items, see \code{\link[sp]{spplot}}.
  }
  %% currently hard-coded for simplicity (rarely needed):
  %% \item{aspect}{the aspect ratio, see \code{\link[lattice]{levelplot}}.
  %%   Defaults to \code{"iso"} for projected coordinates.}
  \item{xlim,ylim}{numeric vectors of length 2 specifying the axis limits.}
}

\value{
  a lattice plot of class
  \code{"\link[lattice:trellis.object]{trellis}"}, but see
  \code{\link[sp]{spplot}}.
}

\author{
  Sebastian Meyer
}

\seealso{
  the central \code{\link{stsplot}}-documentation for an overview of
  plot types, and \code{\link{animate.sts}} for animations of
  \code{"sts"} objects.
}

\examples{
data("measlesWeserEms")

# default plot: total region-specific counts over all weeks
plot(measlesWeserEms, type = observed ~ unit)
stsplot_space(measlesWeserEms)  # the same

# cumulative incidence (per 100'000 inhabitants),
# with region labels and white borders
plot(measlesWeserEms, observed ~ unit,
     population = measlesWeserEms@map$POPULATION / 100000,
     labels = list(labels = "GEN", cex = 0.7, font = 3),
     col = "white", lwd = 2,
     sub = "cumulative incidence (per 100'000 inhabitants)")

# incidence in a particular week, manual color breaks, display total
plot(measlesWeserEms, observed ~ unit, tps = 62,
     population = measlesWeserEms@map$POPULATION / 100000,
     at = c(0, 1, 5),
     total.args = list(x = 0, label = "Overall incidence: "))

# if we had only observed a subset of the regions
plot(measlesWeserEms[,5:11], observed ~ unit,
     gpar.missing = list(col = "gray", lty = 4))
}

\keyword{hplot}
\keyword{spatial}
