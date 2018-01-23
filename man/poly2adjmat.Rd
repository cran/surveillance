\name{poly2adjmat}
\alias{poly2adjmat}
\title{
Derive Adjacency Structure of \code{"SpatialPolygons"}
}

\description{
Wrapping around functionality of the \pkg{spdep} package, this function
computes the symmetric, binary (0/1), adjacency matrix from a
\code{"\linkS4class{SpatialPolygons}"} object.
It essentially applies
\code{\link[spdep]{nb2mat}(\link[spdep]{poly2nb}(SpP, ...), style="B",
  zero.policy=zero.policy)}.
}

\usage{
poly2adjmat(SpP, ..., zero.policy = TRUE)
}

\arguments{
  \item{SpP}{an object inheriting from \code{"\linkS4class{SpatialPolygons}"}.}
  \item{\dots}{arguments passed to \code{\link[spdep]{poly2nb}}.
    Its \code{snap} argument might be particularly useful to handle maps
    with sliver polygons.}
  \item{zero.policy}{logical indicating if islands are allowed, see
    \code{\link[spdep]{nb2mat}}.}
}

\value{
a symmetric numeric indicator matrix of size \code{length(SpP)}^2
representing polygon adjacencies.
}

\author{
(of this wrapper) Sebastian Meyer
}

\seealso{
\code{\link[spdep]{poly2nb}} in package \pkg{spdep}
}

\examples{
if (requireNamespace("spdep")) {
    ## generate adjacency matrix for districts of Bayern and Baden-Wuerttemberg
    data("fluBYBW")
    adjmat <- poly2adjmat(fluBYBW@map)

    ## same as already stored in the neighbourhood slot (in different order)
    stopifnot(all.equal(adjmat,
                        neighbourhood(fluBYBW)[rownames(adjmat),colnames(adjmat)]))

    ## a visual check of the district-specific number of neighbours
    plot(fluBYBW@map)
    text(coordinates(fluBYBW@map), labels=rowSums(adjmat==1), font=2, col=2)

    ## the neighbourhood graph can be plotted with spdep
    plot(spdep::mat2listw(adjmat), coordinates(fluBYBW@map))
}
}

\keyword{spatial}
\keyword{graphs}
