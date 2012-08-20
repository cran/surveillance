\encoding{latin1}
\name{R0}
\alias{R0}
\alias{R0.twinstim}
\alias{R0.simEpidataCS}

\title{Computes basic reproduction numbers from fitted models}

\description{
The S3 generic function \code{R0} defined in package \pkg{surveillance} is intended to
compute basic reproduction numbers from fitted epidemic models.
The package currently defines a method for the \code{"\link{twinstim}"} class, which
computes mean numbers of infections caused by infected individuals depending on the event type
and marks attached to the individual, which contribute to the infection pressure
in the epidemic predictor of that class.
There is also a method for simulated epidemics (internally just a
wrapper for the \code{"twinstim"}-method).
}

\usage{
R0(object, ...)

\method{R0}{twinstim}(object, newevents, trimmed = TRUE, ...)
\method{R0}{simEpidataCS}(object, trimmed = TRUE, ...)
}

\arguments{
  \item{object}{A fitted epidemic model object for which an \code{R0} method exists.}
  \item{newevents}{
    an optional \code{data.frame} of events for which the basic reproduction
    numbers should be calculated.  If omitted, it is calculated for the
    original events from the fit.  In this case, if
    \code{trimmed = TRUE} (the default), the result is just
    \code{object$R0}; however, if \code{trimmed = FALSE}, the model
    environment is required, i.e. \code{object} must have been fitted
    with \code{model = TRUE}.
    
    For the \code{twinstim} method, \code{newevents} must at least
    contain the following columns:
    the \code{time} of the events, the factor variable \code{type}, the
    interaction ranges \code{eps.t} and \code{eps.s}, as well as columns
    for the marks used in the epidemic component of the fitted twinstim
    \code{object} as stored in \code{formula(object)$epidemic}.
    The coding of the variables must of course be the same as used for
    fitting. For \code{trimmed} R0 values, \code{newevents} must additionally
    contain the components \code{.influenceRegion} and, if using the
    \code{Fcircle} trick in the \code{siaf} specification, also
    \code{.bdist} (cf. the hidden columns in the \code{events}
    component of class \code{"epidataCS"}).
  }
  \item{trimmed}{
    logical indicating if the individual reproduction numbers should be
    calculated by integrating the epidemic intensities over the
    observation period and region only (\code{trimmed = TRUE}) or over
    the whole time-space domain R+ x R^2 (\code{trimmed = FALSE}). By
    default, if \code{newevents} is missing, the trimmed \code{R0}
    values stored in \code{object} are returned. Trimming means that
    events near the (spatial or temporal) edges of the observation
    domain have lower reproduction numbers (ceteris paribus) because
    events outside the observation domain are not observed.
  }
  \item{\dots}{additional arguments passed to methods.
    Currently unused for the \code{twinstim} method.}
}

\details{
For the \code{"\link{twinstim}"} class, the individual-specific mean
number \eqn{\mu_j} of infections caused by individual (event) \eqn{j}
inside its theoretical (untrimmed) spatio-temporal range of interaction
given by its \code{eps.t} (\eqn{\epsilon}) and \code{eps.s}
(\eqn{\delta}) values is defined as follows (cf. Meyer et al, 2012):
\deqn{\mu_j = e^{\eta_j} \cdot \int_0^\epsilon g(t) dt \cdot
       \int_{b(\bold{0},\delta)} f(\bold{s}) d\bold{s} .}
Here, \eqn{b(\bold{0},\delta)} denotes the disc centred at (0,0)' with
radius \eqn{\delta}, \eqn{\eta_j} is the epidemic linear predictor,
\eqn{g(t)} is the temporal interaction function, and \eqn{f(\bold{s})}
is the spatial interaction function.
Alternatively, the \code{trimmed} (observed) mean reproduction numbers
are obtain by integrating over the observed infectious domains of the
individuals, i.e. integrate \eqn{f} over the intersection of the
influence region with the observation region \code{W}
(i.e. over \eqn{\{ W \cap b(\bold{s}_j,\delta) \} - \bold{s}_j})
and \eqn{g} over the intersection of the observed infectious period with
the observation period \eqn{(t_0;T]} (i.e. over
\eqn{(0; \min(T-t_j,\epsilon)]}).

(Numerical) Integration is performed exactly as during the fitting of
\code{object}, for instance \code{object$nCub} and
\code{object$nCub.adaptive} are used if necessary.
}

\value{
numeric vector of estimated basic reproduction numbers from the fitted
model \code{object} corresponding to the rows of \code{newevents} (if
supplied) or the original fitted events including events of the prehistory.
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.\cr
  DOI-Link: \url{http://dx.doi.org/10.1111/j.1541-0420.2011.01684.x}
}

\author{Sebastian Meyer}

\examples{
# load the 'imdepi' data and a model fit
data(imdepi)
data(imdepifit)

R0s <- R0(imdepifit)
tapply(R0s, imdepi$events@data[names(R0s), "type"], summary)

R0(imdepifit, newevents=marks(imdepi)[1,], trimmed=FALSE)
}

\keyword{methods}
\keyword{univar}
