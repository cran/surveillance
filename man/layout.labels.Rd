\name{layout.labels}
\alias{layout.labels}
\alias{layout.scalebar}

\title{
Layout Items for \code{spplot}
}

\description{
Generate \code{sp.layout} items for use by \code{\link{spplot}}
or plot these items directly in the traditional graphics system.
Function \code{layout.labels} draws labels at the coordinates of the
spatial object, and \code{layout.scalebar} returns a labeled scale bar.
}

\usage{
layout.labels(obj, labels = TRUE, plot = FALSE)

layout.scalebar(obj, corner = c(0.05, 0.95), scale = 1,
                labels = c(0, scale), height = 0.05,
                pos = 3, ..., plot = FALSE)
}

\arguments{
  \item{obj}{
    an object inheriting from a \code{\linkS4class{Spatial}} class.
  }
  \item{labels}{
    specification of the labels. For \code{layout.labels}:
    \itemize{
      \item a \code{FALSE} or \code{NULL} value omits labels
      (\code{NULL} is returned),
      \item \code{labels = TRUE} uses \code{row.names(obj)},
      \item a character or numeric index for a column of \code{obj@data}
      which contains suitable labels,
      \item a vector of length \code{length(obj)} with labels,
      \item or a list of arguments for \code{\link[lattice]{panel.text}},
      where the optional \code{labels} component follows the same rules
      as above.
    }
    For \code{layout.scalebar}, a character vector of length two giving
    the labels to be put above the left and right ends of the scale bar.
  }
  \item{corner}{
    the location of the scale bar in the unit square, where
    \code{c(0,0)} refers to the bottom left corner. By default, the
    scale bar is placed in the top left corner (with a small buffer).
  }
  \item{scale}{
    the width of the scale bar in the units of \code{\link{proj4string}(obj)}.
    If \code{identical(FALSE, \link{is.projected}(obj))}
    (i.e., \code{obj} has longlat coordinates), \code{scale} is
    interpreted in kilometres.
  }
  \item{height}{
    the height of the scale bar, see \code{\link{layout.scale.bar}}.
  }
  \item{pos}{
    a position specifier for the labels (see \code{\link{text}}).
    By default, the labels are plotted above the scale bar.
  }
  \item{\dots}{
    further arguments for \code{\link[lattice]{panel.text}} (if
    \code{plot = FALSE}) or \code{\link{text}} (if \code{plot = TRUE})
    to change the style of the labels, e.g., \code{cex}, \code{col},
    and \code{font}.
  }
  \item{plot}{
    logical indicating if the layout item should be plotted using the
    traditional graphics system. By default (\code{FALSE}), a list for
    subsequent use by \code{\link{spplot}} is returned.
  }
}

\value{
  For \code{layout.labels}, a single \code{sp.layout} item, which is
  a list with first element \code{"panel.text"} and subsequent elements
  being arguments to that function based on the \code{labels}
  specification.

  For \code{layout.scalebar}, a list of \code{sp.layout} items
  comprising the polygonal scale bar and the labels.

  If these layout functions are called with \code{plot = TRUE},
  the item is plotted directly using traditional graphics functions
  and \code{NULL} is returned.
}

\author{
Sebastian Meyer
}

\examples{
## districts in the Regierungsbezirk Weser-Ems (longlat coordinates)
data("measlesWeserEms")
mapWE <- measlesWeserEms@map
li1 <- layout.labels(mapWE, labels = list(font=2, labels="GEN"))
li2 <- layout.scalebar(mapWE, corner = c(0.05, 0.05), scale = 20,
                       labels = c("0", "20 km"))
spplot(mapWE, zcol = "AREA", sp.layout = c(list(li1), li2),
       col.regions = rev(heat.colors(100)), scales = list(draw = TRUE))

## districts in Bavaria (projected coordinates)
load(system.file("shapes", "districtsD.RData", package = "surveillance"))
bavaria <- districtsD[substr(row.names(districtsD), 1, 2) == "09", ]
sb <- layout.scalebar(bavaria, corner = c(0.75,0.9), scale = 50,
                      labels = c("0", "50 km"), cex = 0.8)
spplot(bavaria, zcol = "POPULATION", sp.layout = sb,
       xlab = "x [km]", ylab = "y [km]", scales = list(draw = TRUE),
       col.regions = rev(heat.colors(100)))

## these layout functions also work in the traditional graphics system
par(mar = c(0,0,0,0))
plot(bavaria, col = "lavender")
layout.scalebar(bavaria, corner = c(0.75, 0.9), scale = 50,
                labels = c("0", "50 km"), plot = TRUE)
layout.labels(bavaria, labels = list(cex = 0.8,
              labels = substr(bavaria$GEN, 1, 3)), plot = TRUE)
}

\keyword{aplot}
\keyword{dplot}
