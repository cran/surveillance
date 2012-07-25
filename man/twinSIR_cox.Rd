\name{twinSIR_cox}
\alias{cox}

\title{
Identify Endemic Components in an Intensity Model
}

\description{
The special function \code{cox} marks terms in formulae of the functions
\code{\link{twinSIR}} and \code{\link{simEpidata}} as endemic components,
i.e. variables acting multiplicatively on the baseline infection intensity.

Technically, this function is implemented as \code{function(x) {x}} and
defined as \dQuote{special} in \code{\link{terms.formula}}.
}

\seealso{
Usage in formulae of functions \code{\link{twinSIR}} and
\code{\link{simEpidata}}.
}

\examples{
# a hypothetic call of the fitting function 'twinSIR' with two epidemic and
# two endemic covariates:
\dontrun{
twinSIR( ~ B1 + B2 + cox(vaccination) + cox(size), data = myEpidata)
}
}

\keyword{internal}
