\name{sumNeighbours}
\alias{sumNeighbours}
\title{Calculates the sum of counts of adjacent areas}
\description{
       Calculates the sum of counts of adjacent units/areas, i.e.
       \eqn{\sum_{j \sim i} y_{j,t}} for all time points \eqn{t} and 
       each unit \eqn{i}, \eqn{t=1,\ldots,n, i=1,\ldots,m}.
}
\usage{
  sumNeighbours(disProgObj)
}
\arguments{
        \item{disProgObj}{Object of class \code{disProg}}
}
\value{matrix of dimension \eqn{n \times m}}
}
\keyword{}


% aggregate.disProg
