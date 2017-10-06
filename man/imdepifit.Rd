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
}
\seealso{
  common methods for \code{"twinstim"} fits,
  exemplified using \code{imdepifit}, e.g.,
  \code{\link{summary.twinstim}}, \code{\link{plot.twinstim}},
  and \code{\link{simulate.twinstim}}
}
\examples{
data("imdepi", "imdepifit")

\dontrun{
## reproduce "imdepifit"
myimdepifit <- twinstim(
    endemic = addSeason2formula(~ offset(log(popdensity)) + I(start/365-3.5),
                                S = 1, period = 365, timevar = "start"),
    epidemic = ~ type + agegrp, siaf = siaf.gaussian(),
    data = imdepi, subset = !is.na(agegrp),
    optim.args = list(control = list(reltol = sqrt(.Machine$double.eps))),
    ## the historical default for reltol is 1e-6, which is rather large
    model = FALSE, cumCIF = FALSE
)
stopifnot(all.equal(imdepifit, myimdepifit))
}
}
\keyword{datasets}
