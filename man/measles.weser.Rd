\encoding{latin1}
\name{measles.weser}
\alias{measles.weser}
\alias{measlesWeserEms}
\docType{data}
\keyword{datasets}

\title{Measles in the Weser-Ems region of Lower Saxony, Germany, 2001-2002}

\description{
  Weekly counts of new measles cases for the 17 administrative
  districts (NUTS-3 level) of the \dQuote{Weser-Ems} region of Lower
  Saxony, Germany, during 2001 and 2002, as reported to the Robert Koch
  institute according to the Infection Protection Act
  (\dQuote{Infektionsschutzgesetz}, \acronym{IfSG}).\cr
  \code{data("measlesWeserEms")} is a corrected version of
  \code{data("measles.weser")} (see Format section below).
  These data are illustrated and analyzed in Meyer et al. (2017, Section 5),
  see \code{vignette("hhh4_spacetime")}.
}

\usage{
data("measles.weser")
data("measlesWeserEms")
}

\format{
  \code{data("measles.weser")} is an object of the old \code{"disProg"}
  class, whereas \code{data("measlesWeserEms")} is of the new class
  \code{"\linkS4class{sts}"}.

  Furthermore, the following updates have been applied for
  \code{data("measlesWeserEms")}:
  \itemize{
    \item it includes the two districts \dQuote{SK Delmenhorst} (03401) and
    \dQuote{SK Wilhemshaven} (03405) with zero counts, which are ignored in
    \code{data("measles.weser")}.

    \item it corrects the time lag error for year 2002 caused by a
    redundant pseudo-week \dQuote{0} with 0 counts only
    (the row \code{measles.weser$observed[53,]} is nonsense).
    
    \item it has one more case attributed to \dQuote{LK Oldenburg}
    (03458) during 2001/W17, i.e., 2 cases instead of 1. This reflects
    the official data as of \dQuote{Jahrbuch 2005}, whereas
    \code{data("measles.weser")} is as of \dQuote{Jahrbuch 2004}.

    \item it contains a map of the region (as a
    \code{"\linkS4class{SpatialPolygonsDataFrame}"}) with the following
    variables:
    \describe{
      \item{\code{GEN}}{district label.}
      \item{\code{AREA}}{district area in m^2.}
      \item{\code{POPULATION}}{number of inhabitants (as of 31/12/2003).}
      \item{\code{vaccdoc.2004}}{proportion with a vaccination card
	among screened abecedarians (2004).}
      \item{\code{vacc1.2004}}{proportion with at least one vaccination
	against measles among abecedarians presenting a vaccination card
	(2004).}
      \item{\code{vacc2.2004}}{proportion of doubly vaccinated
	abecedarians among the ones presenting their vaccination card at
	school entry in the year 2004.}
    }

    \item it uses the correct format for the official district keys,
    i.e., 5 digits (initial 0).

    \item its attached neighbourhood matrix is more general: a distance matrix
    (neighbourhood orders) instead of just an adjacency indicator matrix
    (special case \code{nbOrder == 1}).

    \item population fractions represent data as of 31/12/2003 (\acronym{LSN},
    2004, document \dQuote{A I 2 - hj 2 / 2003}). There are only
    minor differences to the ones used for \code{data("measles.weser")}.
  }
}

\source{
  Measles counts were obtained from the public SurvStat database of the
  Robert Koch institute: \url{https://survstat.rki.de/}.

  A shapefile of Germany's districts as of 01/01/2009 was obtained from
  the German Federal Agency for Cartography and Geodesy
  (\url{https://gdz.bkg.bund.de/}).
  The map of the 17 districts of the \dQuote{Weser-Ems} region
  (\code{measlesWeserEms@map}) is a simplified subset of this
  shapefile using a 30\% reduction via the Douglas-Peucker reduction method
  as implemented at \url{https://MapShaper.org}.

  Population numbers were obtained from the Federal Statistical Office of
  Lower Saxony (\acronym{LSN}):
  \url{https://www.statistik.niedersachsen.de/themenbereiche/bevoelkerung/}

  Vaccination coverage was obtained from the public health department of
  Lower Saxony:
  Nieders\enc{ä}{ae}chsisches Landesgesundheitsamt (2005): Impfreport --
  Durchimpfung von Kindern im Einschulungsalter in Niedersachsen im
  Erhebungsjahrgang 2004. Online available from
  \url{https://www.nlga.niedersachsen.de/gesundheitsberichterstattung/gesundheitsberichte/impfreport/},
  also as an interactive version.
}

\references{
  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
}

\examples{
## old "disProg" object
data("measles.weser")
measles.weser
plot(measles.weser, as.one=FALSE)

## new "sts" object (with corrections)
data("measlesWeserEms")
measlesWeserEms
plot(measlesWeserEms)
}
