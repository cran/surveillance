\name{twinstim_methods}
\alias{print.twinstim}
\alias{summary.twinstim}
\alias{coeflist.twinstim}
\alias{vcov.twinstim}
\alias{logLik.twinstim}
\alias{nobs.twinstim}
\alias{print.summary.twinstim}
\alias{toLatex.summary.twinstim}
\alias{xtable.twinstim}
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
  by the corresponding \code{\link{toLatex}} or \code{\link{xtable}}
  methods.
}

\usage{
\method{print}{twinstim}(x, digits = max(3, getOption("digits") - 3), ...)
\method{summary}{twinstim}(object, test.iaf = FALSE,
        correlation = FALSE, symbolic.cor = FALSE, runtime = FALSE, ...)

\method{coeflist}{twinstim}(x, ...)
\method{vcov}{twinstim}(object, ...)
\method{logLik}{twinstim}(object, ...)
\method{nobs}{twinstim}(object, ...)

\method{print}{summary.twinstim}(x,
      digits = max(3, getOption("digits") - 3), symbolic.cor = x$symbolic.cor,
      signif.stars = getOption("show.signif.stars"), ...)

\method{toLatex}{summary.twinstim}(object,
        digits = max(3, getOption("digits") - 3), eps.Pvalue = 1e-4,
        align = "lrrrr", booktabs = getOption("xtable.booktabs", FALSE),
        withAIC = FALSE, ...)
\method{xtable}{summary.twinstim}(x, caption = NULL, label = NULL,
      align = c("l", "r", "r", "r"), digits = 3,
      display = c("s", "f", "s", "s"), ...,
      ci.level = 0.95, ci.fmt = "\%4.2f", ci.to = "--",
      eps.Pvalue = 1e-4)
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
  \item{runtime}{
    logical. If \code{TRUE}, the summary additionally includes the time
    elapsed and the number of log-likelihood and score function evaluations
    during model fitting.
  }
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance
      stars} are printed for each coefficient.}
  \item{eps.Pvalue}{passed to \code{\link{format.pval}}.}
  \item{booktabs}{logical indicating if the \code{toprule},
  \code{midrule} and \code{bottomrule} commands from the LaTeX package
  \pkg{booktabs} should be used for horizontal lines rather than \code{hline}.}
  \item{withAIC}{logical indicating if the AIC and the log-likelihood of
    the model should be included below the table of coefficients in the
    LaTeX tabular.}
  \item{caption,label,align,display}{see \code{\link{xtable}}.}
  \item{ci.level,ci.fmt,ci.to}{the confidence intervals are calculated
    at level \code{ci.level} and printed using \code{\link{sprintf}}
    with format \code{ci.fmt} and separator \code{ci.to}.}
  \item{\dots}{
    For \code{print.summary.twinstim}, arguments passed to
    \code{\link{printCoefmat}}.\cr
    For all other methods: unused (argument of the generic).
  }
}

\details{
  The estimated coefficients and standard Wald-type confidence intervals
  can be extracted using the default \code{\link{coef}} and
  \code{\link{confint}} methods from package \pkg{stats}.
  Note, however, that there is the useful \code{\link{coeflist}} method to
  list the coefficients by model component.

  The \code{print} and \code{summary} methods allow the compact or comprehensive
  representation of the fitting results, respectively.  The former only prints
  the original function call, the estimated coefficients and the maximum
  log-likelihood value.  The latter prints the whole coefficient matrix
  with standard errors, z- and p-values (see \code{\link{printCoefmat}})
  -- separately for the endemic and the epidemic component -- and 
  additionally the AIC, the achieved log-likelihood, the number of
  log-likelihood and score evaluations, and the runtime.
  They both append a big \dQuote{WARNING}, if the optimization algorithm
  did not converge.

  The \code{toLatex} method is essentially a
  translation of the printed summary table of coefficients to LaTeX
  code (using \pkg{xtable}). However, the \code{xtable} method does a
  different job in that it first converts coefficients to rate ratios
  (RR, i.e., the \code{exp}-transformation) and gives confidence
  intervals for those instead of standard errors and z-values.
  Intercepts and interaction function parameters are ignored by the
  \code{xtable} method.
}

\value{
  The \code{print} methods return their first argument, invisibly, as
  they always should.
  The \code{vcov} method returns the estimated variance-covariance
  matrix of the parameters, which is the inverse of
  \code{object$fisherinfo} (estimate of the \emph{expected} Fisher
  information matrix). This \code{"fisherinfo"} is not always available
  (see \code{\link{twinstim}}), in which case
  \code{object$fisherinfo.observed} is used if available or an error is
  returned otherwise. 
  The \code{logLik} and \code{nobs} methods return the maximum
  log-likelihood value of the model, and the number of events (excluding
  events of the prehistory), respectively.
  
  The \code{summary} method returns a list containing some summary
  statistics of the model, which is nicely printed by the corresponding
  \code{print} method.
  
  The \code{toLatex} method returns a character vector of class
  \code{"Latex"}, each element containing one line of LaTeX code (see
  \code{\link{print.Latex}}).
  The \code{xtable} method returns an object of class
  \code{"\link{xtable}"}. Note that the column name of the confidence
  interval, e.g. \dQuote{95\% CI}, contains the percent symbol that may
  need to be escaped when printing the \code{"xtable"} in the output
  format (see \code{sanitize.text.function} in
  \code{\link{print.xtable}}). This may also hold for row names.
}

\author{
Sebastian Meyer
}

\examples{
# load a fit of the 'imdepi' data, see the example in ?twinstim
data("imdepifit")

# print method
imdepifit

# extract point estimates (in a single vector or listed by model component)
coef(imdepifit)
coeflist(imdepifit)

# variance-covariance matrix of endemic parameters
# (inverse of expected Fisher information)
unname(vcov(imdepifit)[1:4,1:4])

# the default confint() method may be used for Wald CI's
confint(imdepifit, parm="e.typeC", level=0.95)

# log-likelihood and AIC of the fitted model
logLik(imdepifit)
AIC(imdepifit)
nobs(imdepifit)

# produce a summary with parameter correlations and runtime information
(s <- summary(imdepifit, correlation=TRUE, symbolic.cor=TRUE, runtime=TRUE))

# create LaTeX code of coefficient table
toLatex(s, withAIC=FALSE)

# or using the xtable-method (which produces rate ratios)
xtable(s)
}

\keyword{methods}
\keyword{print}
\keyword{htest}
