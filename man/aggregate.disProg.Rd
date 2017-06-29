\name{aggregate.disProg}
\alias{aggregate.disProg}
\title{Aggregate the observed counts}
\description{
       Aggregates the observed counts for a multivariate
       \code{disProgObj} over the units. Future versions 
       of \code{surveillance} will also allow for time aggregations etc. 
}
\usage{
  \method{aggregate}{disProg}(x,\dots)
}
\arguments{
        \item{x}{Object of class \code{disProg}}
        \item{\dots}{not used at the moment}
}
\value{\item{x}{univariate \code{disProg} object with aggregated 
     counts and respective states for each time point.}
}
\keyword{hplot}
\examples{
data(ha)
plot(aggregate(ha))
}
