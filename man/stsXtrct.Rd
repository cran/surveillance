\name{[-methods}
\docType{methods}
\alias{[-methods}
\title{Methods for "[": Extraction or Subsetting in Package 'surveillance'}
%
% This is stolen from the Matrix package. 
%\alias{[,sts,missing,missing,ANY-method}
%\alias{[,sts,missing,index,missing-method}
%\alias{[,sts,index,missing,missing-method}
\alias{[,sts-method}
%
\description{
  Methods for \code{"["}, i.e., extraction or subsetting of the
  \code{sts} class in package \pkg{surveillance}.

  Note that \code{[<--methods} methods (i.e. subassigments) are
  currently not supported.

  \code{drop} is always \code{FALSE}.
}
\section{Methods}{
  There are more than these:
  \describe{
    \item{x = "sts", i = "missing", j = "missing", drop= "ANY"}{ ... }
    \item{x = "sts", i = "numeric", j = "missing", drop= "missing"}{ ... }
    \item{x = "sts", i = "missing", j = "numeric", drop= "missing"}{ ... }
  }
}

\examples{
data(ha)
#Convert to S4 object
shp <- system.file("shapes/berlin.shp",package="surveillance")
has4 <- aggregate(disProg2sts(ha,map=readShapePoly(shp,IDvar="SNAME")),nfreq=13)

#A suite of of simple tests (inspired by the Matrix package)
stopifnot(identical(has4, has4[]))


plot(has4[, 3])       # Single series
plot(has4[1:30, 3])   # Somewhat shorter

#Counts at time 20
plot(has4[20, ],type = observed ~1 |unit)

}
\keyword{methods}
\keyword{array}
