\name{[,sts-methods}
\docType{methods}
\title{Extraction and Subsetting of \pkg{sts} objects}
%
% This is stolen from the Matrix package. 
%\alias{[,sts,missing,missing,ANY-method}
%\alias{[,sts,missing,index,missing-method}
%\alias{[,sts,index,missing,missing-method}
\alias{[,sts-method}
\alias{[,sts,ANY,ANY,ANY-method}
%
\description{
  Methods for \code{"["}, i.e., extraction or subsetting of the
  \code{"\linkS4class{sts}"} class in package \pkg{surveillance}.

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
data("ha.sts")
haagg <- aggregate(ha.sts, nfreq=13)

#A suite of of simple tests (inspired by the Matrix package)
#stopifnot(identical(ha.sts, ha.sts[]))


plot(haagg[, 3])       # Single series
plot(haagg[1:30, 3])   # Somewhat shorter

#Counts at time 20
plot(haagg[20, ], type = observed ~1 |unit)

}
\keyword{methods}
\keyword{array}
