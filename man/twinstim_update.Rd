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
\method{update}{twinstim}(object, endemic, epidemic, optim.args, model,
       ..., evaluate = TRUE)
}

\arguments{
  \item{object}{a previous \code{"twinstim"} fit.}
  \item{endemic, epidemic}{changes to the formulae -- see
    \code{\link{update.formula}} and \code{\link{twinstim}}.}
  \item{optim.args}{see \code{\link{twinstim}}. If a list, it will
    modify the original \code{optim.args} using \code{\link{modifyList}}.}
  \item{model}{see \code{\link{twinstim}}. If this is the only argument
    to update, re-fitting is cleverly circumvented.}
  \item{\dots}{Additional arguments to the call, or arguments with changed
    values.}
  \item{evaluate}{If \code{TRUE} (default), evaluate the new call else
    return the call.}
}

\value{
  If \code{evaluate = TRUE} the re-fitted object, otherwise the updated call.
}

\author{
  Sebastian Meyer
  
  Inspiration and some pieces of code originate from
  \code{update.default} by the R Core Team.
}

\seealso{
\code{\link{update.default}}
}

\examples{
if (surveillance.options("allExamples"))
{
  data(imdepi)
  data(imdepifit)

  ## enrich the fit by the model environment
  imdepifit <- update(imdepifit, model=TRUE)
  ## this enables, e.g., intensityplot()s and untrimmed R0 estimates

  ## add another epidemic covariate
  imdepifit2 <- update(imdepifit, epidemic = ~. + log(popdensity),
    # must also update vector of start values!
    optim.args=list(par=c(head(coef(imdepifit),-1),0,tail(coef(imdepifit),1)),
                    fixed="e.siaf.1"), # for reasons of speed
    cumCIF=FALSE)
  
  ## compare with AIC
  AIC(imdepifit, imdepifit2)
}
}

\keyword{models}
\keyword{methods}
