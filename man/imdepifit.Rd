\name{imdepifit}
\alias{imdepifit}
\docType{data}
\title{
Example \code{twinstim} Fit for the \code{imdepi} Data
}
\description{
  \code{data("imdepifit")} is a \code{\link{twinstim}} model
  fitted to the \code{\link{imdepi}} data.
}
\usage{data("imdepifit")}
\format{
  an object of class \code{"\link{twinstim}"}
  obtained from the following call using \code{data(imdepi)}:
  
  \Sexpr[stage=build,results=verbatim]{
    data(imdepifit, package="surveillance")
    writeLines(deparse(imdepifit$call, width.cutoff = 25))
  }
}
\seealso{
  common methods for \code{"twinstim"} fits,
  exemplified using \code{imdepifit}, e.g.,
  \code{\link{summary.twinstim}}, \code{\link{plot.twinstim}},
  and \code{\link{simulate.twinstim}}
}
\examples{
data("imdepi", "imdepifit")

## how this fit was obtained
imdepifit$call
\dontshow{if (surveillance.options("allExamples"))
## reproduce "imdepifit"
stopifnot(all.equal(imdepifit, eval(imdepifit$call)))
}
}
\keyword{datasets}
