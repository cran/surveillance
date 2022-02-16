\name{aggregate.disProg}
\alias{aggregate.disProg}
\title{Aggregate a \code{disProg} Object}
\description{
       Aggregates the observed counts of a multivariate
       \code{disProg} object over the units.
}
\usage{
  \method{aggregate}{disProg}(x,\dots)
}
\arguments{
        \item{x}{Object of class \code{disProg}}
        \item{\dots}{not used at the moment}
}
\value{univariate \code{disProg} object with aggregated 
     counts and respective states for each time point.}
\keyword{internal} % disProg is soft-deprecated
\examples{
data(ha)
dim(ha$observed)
dim(aggregate(ha)$observed)
}
