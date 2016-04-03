\name{[,sts-methods}
\docType{methods}
\title{Extraction and Subsetting of \code{"sts"} Objects}
\alias{[,sts-method}
\alias{[,sts,ANY,ANY,ANY-method}
\description{
  \code{"["}-methods, i.e., extraction or subsetting of the
  \code{"\linkS4class{sts}"} class in package \pkg{surveillance}.

  Note that \code{[<-} methods (i.e. subassignments) are
  currently not supported.

  The \code{drop} argument is always \code{FALSE}.
}

\examples{
data("ha.sts")
haagg <- aggregate(ha.sts, nfreq=13)

plot(haagg[, 3])       # Single series
plot(haagg[1:30, 3])   # Somewhat shorter

#Counts at time 20
plot(haagg[20, ], type = observed ~ unit)
}

\keyword{methods}
\keyword{array}
