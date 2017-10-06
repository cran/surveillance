\name{twinstim_update}
\alias{update.twinstim}
\title{
\code{update}-method for \code{"twinstim"}
}

\description{
Update and (by default) re-fit a \code{"twinstim"}. This method is
especially useful if one wants to add the \code{model} environment
(which is required for some methods) to a fitted model object a posteriori.
}

\usage{
\method{update}{twinstim}(object, endemic, epidemic,
       control.siaf, optim.args, model,
       ..., use.estimates = TRUE, evaluate = TRUE)
}

\arguments{
  \item{object}{a previous \code{"twinstim"} fit.}
  \item{endemic, epidemic}{changes to the formulae -- see
    \code{\link{update.formula}} and \code{\link{twinstim}}.}
  \item{control.siaf}{a list (see \code{\link{twinstim}}) to replace
    the given elements in the original \code{control.siaf} list.
    If \code{NULL}, the original list of control arguments is removed
    from the call, i.e., the defaults are used in \code{twinstim}.}
  \item{optim.args}{see \code{\link{twinstim}}. If a list, it will
    modify the original \code{optim.args} using \code{\link{modifyList}}.}
  \item{model}{see \code{\link{twinstim}}. If this is the only argument
    to update, re-fitting is cleverly circumvented.
    Enriching the fit by the model environment is, e.g., required for
    \code{\link{intensityplot.twinstim}}.}
  \item{\dots}{Additional arguments to the call, or arguments with changed
    values.\cr
    If \code{start} values are specified, they need to be in the same
    format as in the original call \code{object$call$start},
    which is either a named list of named numeric vectors or a named
    numeric vector; see the argument description in \code{\link{twinstim}}.}
  \item{use.estimates}{logical indicating if the estimates of
  \code{object} should be used as initial values for the new fit
  (in the \code{start} argument of \code{twinstim}). Defaults to \code{TRUE}.}
  \item{evaluate}{If \code{TRUE} (default), evaluate the new call else
    return the call.}
}

\value{
  If \code{evaluate = TRUE} the re-fitted object, otherwise the updated call.
}

\author{
  Sebastian Meyer
  
  Inspiration and some pieces of code originate from
  \code{\link{update.default}} by the R Core Team.
}

\seealso{
\code{\link{update.default}}
}

\examples{
data("imdepi", "imdepifit")

## add another epidemic covariate
## (but fix siaf-parameter so that this example runs quickly)
imdepifit2 <- update(imdepifit, epidemic = ~. + log(popdensity),
                     optim.args = list(fixed="e.siaf.1"))

## compare by AIC
AIC(imdepifit, imdepifit2)
}

\keyword{models}
\keyword{methods}
