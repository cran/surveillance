\name{disProg2sts}
\alias{disProg2sts}
\title{Convert disProg object to sts}
\description{
  A small helper function to convert a \code{disProg} object to become an
  object of the S4 class \code{sts}. In the future the \code{sts} should
  replace the \code{disProg} class.
}
\usage{
   disProg2sts(disProgObj, map=NULL)
}
\arguments{
  \item{disProgObj}{object of class disProg}
  \item{map}{\code{SpatialPolygonsDataFrame} object containing the
    map visualization}
}
\value{
  an object of class \code{sts}
}
\seealso{
  \code{\link{sts-class}}
}
\examples{
  data(ha)
  print(disProg2sts(ha))
}
\keyword{}

