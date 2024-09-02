\name{momo}
\alias{momo}
\docType{data}
\encoding{latin1}

\title{Danish 1994-2008 all-cause mortality data for eight age groups}
\description{
  Weekly number of all cause mortality from 1994-2008 in each
  of the eight age groups <1, 1-4, 5-14, 15-44, 45-64, 65-74,
  75-84 and 85+ years, see \enc{Höhle}{Hoehle} and Mazick (2010).
}
\usage{data(momo)}
\format{
  An object of class \code{"\linkS4class{sts}"} containing the weekly
  number of all-cause deaths in Denmark, 1994-2008 (782 weeks), for each
  of the eight age groups <1, 1-4, 5-14, 15-44, 45-64, 65-74,
  75-84 and 85+ years. A special feature of the EuroMOMO data is that
  weeks follow the ISO 8601 standard, which can be
  handled by the \code{"sts"} class.

  The \code{population} slot of the \code{momo} object contains the
  population size in each of the eight age groups.
  These are yearly data obtained from the StatBank Denmark.
}
\source{
  \emph{European monitoring of excess mortality for public health action}
  (EuroMOMO) project. \url{https://www.euromomo.eu/}.
  
  Department of Epidemiology, Statens Serum Institute, Copenhagen, Denmark
  StatBank Denmark, Statistics Denmark, \url{https://www.statistikbanken.dk/}
}
\examples{
data("momo")
momo

## show the period 2000-2008 with customized x-axis annotation
## (this is Figure 1 in Hoehle and Mazick, 2010)
oopts <- surveillance.options("stsTickFactors" = c("\%G" = 1.5, "\%Q"=.75))
plot(momo[year(momo) >= 2000,], ylab = "", xlab = "Time (weeks)",
     par.list = list(las = 1), col = c(gray(0.5), NA, NA),
     xaxis.tickFreq = list("\%G"=atChange, "\%Q"=atChange),
     xaxis.labelFreq = list("\%G"=atChange), xaxis.labelFormat = "\%G")
surveillance.options(oopts)


if (surveillance.options("allExamples")) {

## stratified monitoring from 2007-W40 using the Farrington algorithm
phase2 <- which(epoch(momo) >= "2007-10-01")
momo2 <- farrington(momo, control = list(range=phase2, alpha=0.01, b=5, w=4))
print(colSums(alarms(momo2)))
plot(momo2, col = c(8, NA, 4), same.scale = FALSE)

## stripchart of alarms (Figure 5 in Hoehle and Mazick, 2010)
plot(momo2, type = alarm ~ time, xlab = "Time (weeks)", main = "",
     alarm.symbol = list(pch=3, col=1, cex=1.5))

}
}
\references{
  \enc{Höhle}{Hoehle}, M. and Mazick, A. (2010). Aberration detection in R
  illustrated by Danish mortality monitoring. In T. Kass-Hout and X.
  Zhang (eds.), \emph{Biosurveillance: A Health Protection Priority},
  chapter 12. Chapman & Hall/CRC.\cr
  Preprint available at
  \url{https://staff.math.su.se/hoehle/pubs/hoehle_mazick2009-preprint.pdf}
}
\keyword{datasets}
