\name{disProg2sts}
\alias{disProg2sts}
\alias{sts2disProg}
\title{Convert disProg object to sts and vice versa}
\description{
  A small helper function to convert a \code{disProg} object to become
  an object of the S4 class \code{sts} and vice versa. In the future the
  \code{sts} should replace the \code{disProg} class, but for now this
  function allows for conversion between the two formats.
}
\usage{
   disProg2sts(disProgObj, map=NULL)
   sts2disProg(sts)
}
\arguments{
  \item{disProgObj}{an object of class \code{"disProg"}}
  \item{map}{an optional \code{"SpatialPolygons"} object}
  \item{sts}{an object of class \code{"sts"} to convert}
}
\value{
  an object of class \code{"sts"} or \code{"disProg"}, respectively.
}
\seealso{
  \code{\link{sts-class}}
}
\examples{
  data(ha)
  print(disProg2sts(ha))
  class(sts2disProg(disProg2sts(ha)))
}
\keyword{utilities}
