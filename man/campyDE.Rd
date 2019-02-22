\encoding{latin1}
\name{campyDE}
\alias{campyDE}
\docType{data}
\title{Campylobacteriosis and Absolute Humidity in Germany 2002-2011}
\description{
  Weekly number of reported campylobacteriosis cases in Germany,
  2002-2011, together with the corresponding absolute humidity (in g/m^3)
  that week. The absolute humidity was computed according to the
  procedure by Dengler (1997) using the means of representative weather
  station data from the German Climate service.
}
\usage{
data(campyDE)
}
\format{
  A \code{data.frame} containing the following columns
  \describe{
    \item{\code{date}}{\code{Date} instance containing the Monday of the
      reporting week.}
    \item{\code{case}}{Number of reported cases that week.}
    \item{\code{state}}{Boolean indicating whether there is external knowledge
      about an outbreak that week}
    \item{\code{hum}}{Mean absolute humidity (in g/m^3) of that week as
      measured by a single representative weather station.}
    \item{\code{l1.hum}-\code{l5.hum}}{Lagged version (lagged by 1-5) of
      the \code{hum} covariate.}
    \item{newyears}{Boolean indicating whether the reporting week
      corresponds to the first two weeks of the year (TRUE) or not
      (FALSE). Note: The first week of a year is here defined as the first
      reporting week, which has its corresponding Monday within new year.}
    \item{christmas}{Boolean indicating whether the reporting week
      corresponds to the last two weeks of the year (TRUE) or not
      (FALSE). Note: This are the first two weeks before the
      \code{newyears} weeks.} 
    \item{O104period}{Boolean indicating whether the reporting week
      corresponds to the W21-W30 period of increased gastroenteritis
      awareness during the O104:H4 STEC outbreak.}
  }
}
\source{
  The data on campylobacteriosis cases have been queried from the
  Survstat@RKI database of the German Robert Koch Institute
  (\url{https://survstat.rki.de/}).

  Data for the computation of absolute humidity were obtained from the
  German Climate Service (Deutscher Wetterdienst), Climate data of
  Germany, available at \url{http://www.dwd.de}.

  A complete data description and an analysis of the data can be found
  in Manitz and \enc{Höhle}{Hoehle} (2013).
}
\references{
  Manitz, J. and \enc{Höhle}{Hoehle}, M. (2013):
  Bayesian outbreak detection algorithm for monitoring reported cases of
  campylobacteriosis in Germany.
  Biometrical Journal, 55(4), 509-526.
}
\examples{
# Load the data
data("campyDE")

# O104 period is W21-W30 in 2011
stopifnot(all(campyDE$O104period == (
  (campyDE$date >= as.Date("2011-05-23")) &
  (campyDE$date < as.Date("2011-07-31"))
)))

# Make an sts object from the data.frame
cam.sts <- sts(epoch=campyDE$date, observed=campyDE$case, state=campyDE$state)

# Plot the result
plot(cam.sts)
}
\keyword{datasets}
