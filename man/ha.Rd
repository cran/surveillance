\name{ha}
\alias{ha}
\alias{ha.sts}
\docType{data}
\title{Hepatitis A in Berlin}

\description{
  Number of Hepatitis A cases among adult male (age>18) in
  Berlin 2001-2006. An increase is seen during 2006
}

\usage{
data("ha")
data("ha.sts")
}

\format{
  \code{ha} is a \code{disProg} object containing \eqn{290\times 12}{290 x 12}
  observations starting from week 1 in 2001 to week 30 in 2006.
  \code{ha.sts} is generated from \code{ha} by the converter function
  \code{\link{disProg2sts}} using a shape file of Berlin (see Examples).
}

\source{
Robert Koch-Institut: SurvStat: \url{https://survstat.rki.de/};
Queried on 25 August 2006.

Robert Koch Institut, Epidemiologisches Bulletin 33/2006, p.290.
}

\examples{
## deprecated "disProg" object
data("ha")
ha
plot(aggregate(ha))

## new-style "sts" object
data("ha.sts")
ha.sts
plot(ha.sts, type = observed ~ unit, labels = TRUE)

## conversion of the old "disProg" object 'ha' to the new S4 class "sts"
\dontrun{
shpfile <- system.file("shapes/berlin.shp", package="surveillance")
ha.sts <- disProg2sts(ha, map = maptools::readShapePoly(shpfile,IDvar="SNAME"))
## in data("ha.sts"), German umlauts in 'ha.sts@map@data$BEZIRK'
## have been replaced for compatibility
}
}

\keyword{datasets}
