\name{twinstim_simEndemicEvents}
\alias{simEndemicEvents}
\title{
  Quick Simulation from an Endemic-Only \code{twinstim}
}
\description{
In \emph{endemic-only} \code{\link{twinstim}} models, the conditional
intensity is a piecewise constant function independent from the history
of the process. This allows for a much more efficient simulation
algorithm than via Ogata's modified thinning as in the general
\code{\link{simulate.twinstim}} method.
}
\usage{
simEndemicEvents(object, tiles)
}
\arguments{
  \item{object}{
    an object of class \code{"\link{twinstim}"} (with the \code{model}
    component retained; otherwise try
    \code{object <- \link[=update.twinstim]{update}(object, model = TRUE)}).
  }
  \item{tiles}{
    an object inheriting from \code{"\linkS4class{SpatialPolygons}"},
    which represents the tiles of the original data's \code{stgrid}
    (see, e.g., \code{levels(environment(object)$gridTiles)}).
  }
}
\value{
  a \code{\linkS4class{SpatialPointsDataFrame}}
}
\author{
  Sebastian Meyer
}
\seealso{
  the general simulation method \code{\link{simulate.twinstim}}
}
\examples{
data("imdepi", "imdepifit")
load(system.file("shapes", "districtsD.RData", package="surveillance"))

## Fit an endemic-only twinstim()
m_noepi <- update(imdepifit, epidemic = ~0, siaf = NULL, model = TRUE)

## Simulate events from the above endemic model
set.seed(1)
s1 <- simEndemicEvents(m_noepi, tiles = districtsD)
class(s1)  # just a "SpatialPointsDataFrame"
summary(s1)
plot(s1, col = s1$type, cex = 0.5); plot(imdepi$W, lwd = 2, add = TRUE)

\dontrun{
## the general simulation method takes several seconds
s0 <- simulate(m_noepi, seed = 1, data = imdepi, tiles = districtsD)
class(s0)  # gives a full "simEpidataCS" with several methods applicable
methods(class = "epidataCS")
plot(s0, "time")
plot(s0, "space", points.args = list(pch = 3), lwd = 2)
}
}
\keyword{datagen}
\keyword{models}
