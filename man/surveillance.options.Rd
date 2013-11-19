\name{surveillance.options}
\alias{surveillance.options}
\alias{reset.surveillance.options}

\title{Options of the \pkg{surveillance} Package}

\description{
  Query, set or reset options specific to the \pkg{surveillance}
  package, similar to what \code{\link{options}} does for global settings.
}

\usage{
surveillance.options(...)
reset.surveillance.options()
}

\arguments{
  \item{\dots}{
    Either empty, or a sequence of option names (as strings),
    or a sequence of \code{name=value} pairs, or a named list of
    options. Available options are:
    \describe{
      \item{gpclib:}{
        Logical flag indicating whether \pkg{gpclib}, the General
	Polygon Clipping Library for \R, which has a restricted license
	(commercial use prohibited), 
	may be used. This is no longer required since package
	\pkg{surveillance} has switched to alternatives such as
	\pkg{polyclip} and \pkg{rgeos} for generating \code{"epidataCS"}
	objects by \code{as.epidataCS} or \code{simEpidataCS}.
	However, for \code{\link{unionSpatialPolygons}} and
	\code{\link{intersectPolyCircle.gpc.poly}}, using \pkg{gpclib}
	is still an option (mainly for backwards compatibility).
	The default setting is \code{FALSE}.
      }

      \item{allExamples:}{
	Logical flag mainly for CRAN-compatibility, i.e. to prevent
	cumbersome computations in help file examples from being run by
	CRAN servers. This option defaults to \code{TRUE} unless the
	environment variable \env{_R_CHECK_TIMINGS_} is set when
	attaching the \pkg{surveillance} package.
      }
    }
  }
}

\value{
  \code{reset.surveillance.options} reverts all options to their default
  values and (invisibly) returns these in a list.

  For \code{surveillance.options}, the following holds:
  \itemize{
    \item If no arguments are given, the current values of all package options are
    returned in a list.
    \item If one option name is given, the current value of this
    option is returned (\emph{not} in a list, just the value).
    \item If several option names are given, the current values of these options
    are returned in a list.
    \item If \code{name=value} pairs are given, the named options
    are set to the given values, and the \emph{previous} values of
    these options are returned in a list.
  }
}

\author{
  Sebastian Meyer,
  inspired by the implementation of \code{spatstat.options()} in the \pkg{spatstat}
  package by Adrian Baddeley and Rolf Turner.
}

\examples{
surveillance.options()
}

\keyword{environment}
