\name{layout.labels}
\alias{layout.labels}
\title{
Layout Item for Feature Labels in \code{spplot}
}

\description{
Generate an \code{sp.layout} item for use in \code{\link{spplot}} in
order to draw labels at the coordinates of the spatial object.
}

\usage{
layout.labels(obj, labels = TRUE)
}

\arguments{
  \item{obj}{
    an object inheriting from a \code{\linkS4class{Spatial}} class.
  }
  \item{labels}{
    specification of the labels:
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
  }
}

\value{
  an \code{sp.layout} item for \code{\link{spplot}}, which is a list.
  Its first element is \code{"panel.text"} and subsequent elements are
  arguments to that function based on the \code{labels} specification.
}

\author{
Sebastian Meyer
}

\examples{
data("measlesWeserEms")
(li <- layout.labels(measlesWeserEms@map, labels = list(font=2, labels="GEN")))

## example usage in spplot()
spplot(measlesWeserEms@map, zcol = "AREA", sp.layout = list(li),
       col.regions = rev(heat.colors(100)))
}

\keyword{aplot}
\keyword{dplot}
