\encoding{latin1}
\name{twinstim_iaf}
\alias{siaf.constant}
\alias{siaf.gaussian}
\alias{siaf.lomax}
\alias{tiaf.constant}
\alias{tiaf.exponential}

\title{
  Temporal and spatial interaction functions for \code{twinstim}
}
\description{
  Predefined functional forms for the spatial and temporal dispersal
  (\eqn{f} and \eqn{g}, respectively) for use in \code{twinstim} models as
  described in Meyer et al. (2012). 
}
\usage{
siaf.constant()
siaf.gaussian(nTypes, logsd = TRUE, density = FALSE, effRangeMult = 6,
              validpars = NULL)
siaf.lomax(nTypes = 1, logpars = TRUE, density = FALSE, effRangeProb = 0.99,
           validpars = NULL)

tiaf.constant()
tiaf.exponential(nTypes)
}
\arguments{
  \item{nTypes}{Number of different types, i.e. \eqn{|\mathcal{K}|}{|K|}}
  \item{logsd}{Parametrize standard deviation parameter
  \eqn{\sigma>0}{sigma>0} on log scale to enforce parameter
  constraint. Default: \code{TRUE}}
  \item{density}{logical. If \code{TRUE}, a density is used as the
    interaction function: for \code{siaf.gaussian}, the density of the
    bivariate, isotropic normal distribution, for \code{siaf.lomax},
    \eqn{f(||\bold{s}||)} is the density of the (univariate) Lomax
    distribution (\eqn{f(\bold{s})} is not a bivariate density here).
    However, in this case, the magnitude of the
    interaction function depends on the parameters of the distribution,
    which might interfere with the identification of an epidemic
    intercept. Hence, by default (\code{density = FALSE}), the density
    is normalized to unity at the origin, i.e.
    \eqn{f(\bold{s}) / f(\bold{0})} is used.}
  \item{effRangeMult}{Determines the effective range for numerical integration
  in terms of multiples of \eqn{\sigma}, i.e. with \code{effRangeMult=6} numerical
  integration only considers the 6-sigma area around the event instead of the whole
  observation region W.}
  \item{validpars}{See \code{\link{twinstim}}. If \code{logsd = FALSE}, you should
  either use constrained optimisation (L-BFGS-B) or set \code{validpars} to
  \code{function (pars) pars > 0}.}
  \item{logpars}{logical. Parametrize the Lomax distribution with log(scale) and
    log(shape) to enforce positive parameters. The alternative
    (\code{FALSE}) is not implemented.}
  \item{effRangeProb}{Determines the effective range for numerical integration
  in terms of the corresponding quantile of the univariate Lomax distribution.}
}
\value{
A list containing the interaction function, its derivative, ...
%TODO
}
\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}

  Meyer, S. (2010):
  Spatio-Temporal Infectious Disease Epidemiology based on Point Processes.
  Master's Thesis, Ludwig-Maximilians-Universit\enc{ä}{ae}t
  M\enc{ü}{ue}nchen.\cr
  Available as \url{http://epub.ub.uni-muenchen.de/11703/}
}
\author{
  S. Meyer with documentation contributions by M. \enc{Höhle}{Hoehle}.
}
\seealso{
\code{\link{twinstim}}.
}
\examples{
#Finetype indenepdent spatial interaction function used in the twinstim paper
siaf_log1 <- siaf.gaussian(1, logsd = TRUE, density = FALSE, effRangeMult = 6)
}
\keyword{models}
\keyword{utilities}

