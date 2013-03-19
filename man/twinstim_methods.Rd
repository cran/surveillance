\encoding{latin1}
\name{twinstim_methods}
\alias{print.twinstim}
\alias{summary.twinstim}
\alias{vcov.twinstim}
\alias{logLik.twinstim}
\alias{nobs.twinstim}
\alias{print.summary.twinstim}
\alias{toLatex.summary.twinstim}
\alias{xtable.summary.twinstim}

\title{
Print, Summary and Extraction Methods for \code{"twinstim"} Objects
}

\description{
  Besides \code{\link{print}} and \code{\link{summary}} methods there
  are also some standard extraction methods defined for objects of class
  \code{"twinstim"}: \code{\link{vcov}}, \code{\link{logLik}}, and
  \code{\link{nobs}}. This
  also enables the use of, e.g., \code{\link{confint}} and
  \code{\link{AIC}}. The model \code{summary} can be exported to LaTeX
  by the corresponding \code{\link{toLatex}} method.
}

\usage{
\method{print}{twinstim}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{twinstim}(object, test.iaf = FALSE,
        correlation = FALSE, symbolic.cor = FALSE, ...)

\method{vcov}{twinstim}(object, ...)
\method{logLik}{twinstim}(object, ...)
\method{nobs}{twinstim}(object, ...)

\method{print}{summary.twinstim}(x,
      digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), ...)

\method{toLatex}{summary.twinstim}(object,
        digits = max(3, getOption("digits") - 3), eps.Pvalue = 1e-4,
        align = "lrrrr", withAIC = FALSE, ...)
\method{xtable}{summary.twinstim}(x, caption = NULL, label = NULL,
      align = c("l", "r", "r", "r"), digits = 3,
      display = c("s", "f", "s", "s"),
      ci.level = 0.95, ci.fmt = "\%4.2f", ci.to = "--",
      eps.pvalue = 1e-4, ...)
}

\arguments{
  \item{x, object}{an object of class \code{"twinstim"} or
    \code{"summary.twinstim"}, respectively.}
  \item{digits}{
    integer, used for number formatting with \code{signif()}.  Minimum number of
    significant digits to be printed in values.
  }
  \item{test.iaf}{logical indicating if the simple Wald z- and p-values
    should be calculated for parameters of the interaction functions
    \code{siaf} and \code{tiaf}.
    Because it is often invalid or meaningless to do so, the default is
    \code{FALSE}.
  }
  \item{correlation}{
    logical. If \code{TRUE}, the correlation matrix of the estimated parameters
    is returned and printed.
  }
  \item{symbolic.cor}{
    logical. If \code{TRUE}, print the correlations in a symbolic form (see
    \code{symnum}) rather than as numbers.
  }
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance
      stars} are printed for each coefficient.}
  \item{eps.Pvalue}{see \code{\link{printCoefmat}}.}
  \item{withAIC}{logical indicating if the AIC and the log-likelihood of
    the model should be included below the table of coefficients in the
    LaTeX tabular.}
  \item{caption,label,align,display}{see \code{\link{xtable}}.}
  \item{ci.level,ci.fmt,ci.to}{the confidence intervals are calculated
    at level \code{ci.level} and printed using \code{\link{sprintf}}
    with format \code{ci.fmt} and seperator \code{ci.to}.}
  \item{eps.pvalue}{argument passed to \code{\link{formatPval}}.}
  \item{\dots}{
    For \code{print.summary.twinstim}, arguments passed to
    \code{\link{printCoefmat}}.\cr
    For all other methods: unused (argument of the generic).
  }
}

\details{
  The \code{print} and \code{summary} methods allow the compact or comprehensive
  representation of the fitting results, respectively.  The former only prints
  the original function call, the estimated coefficients and the maximum
  log-likelihood value.  The latter prints the whole coefficient matrix
  with standard errors, z- and p-values (see \code{\link{printCoefmat}})
  -- separately for the endemic and the epidemic component -- and 
  additionally the AIC, the number of log-likelihood and score
  evaluations, and the runtime.  They both append a big
  \dQuote{WARNING}, if the optimization algorithm did not converge.

  The estimated coefficients may be extracted by using the default
  \code{coef}-method from package \pkg{stats}.
}

\value{
  The \code{print} methods return their first argument, invisibly, as
  they always should. The \code{vcov}, \code{logLik}, and \code{nobs}
  methods return the estimated variance-covariance
  matrix of the parameters (here, the inverse of the estimate of the
  expected Fisher information matrix), the maximum log-likelihood
  value of the model, and the number of events (excluding events of the
  pre-history), respectively.
  The \code{summary} method returns a list containing some summary
  statistics of the model, which is nicely printed by the corresponding
  \code{print} method.
  The \code{toLatex} method returns a character vector of class
  \code{"Latex"}, each element containing one line of LaTeX code (see
  \code{\link{print.Latex}}). 
}

\author{
Sebastian Meyer
}

\examples{
# load a fit of the 'imdepi' data, see the example in ?twinstim
data(imdepifit)

# print method
imdepifit

# extract point estimates
coef(imdepifit)

# variance-covariance matrix of endemic parameters
# (inverse of expected Fisher information)
unname(vcov(imdepifit)[1:4,1:4])

# the default confint() method may be used for Wald CI's
confint(imdepifit, parm="e.typeC", level=0.95)

# log-likelihood and AIC of the fitted model
logLik(imdepifit)
AIC(imdepifit)
nobs(imdepifit)

# summary method
summary(imdepifit, test.iaf=FALSE, correlation=TRUE, symbolic.cor=TRUE)

# create LaTeX code of coefficient table
toLatex(summary(imdepifit), withAIC=FALSE)
}

\keyword{methods}
\keyword{print}
\keyword{htest}
