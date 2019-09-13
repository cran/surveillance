\encoding{latin1}
\docType{data}
\name{imdepi}
\alias{imdepi}

\title{
  Occurrence of Invasive Meningococcal Disease in Germany
}

\description{
  \code{imdepi} contains data on the spatio-temporal location of 636
  cases of invasive meningococcal disease (IMD) caused by the two most
  common meningococcal finetypes in Germany, \samp{B:P1.7-2,4:F1-5} (of
  serogroup B) and \samp{C:P1.5,2:F3-3} (of serogroup C).
}

\usage{
data("imdepi")
}

\format{
  \code{imdepi} is an object of class
  \code{"\link{epidataCS}"} (a list with components \code{events},
  \code{stgrid}, \code{W} and \code{qmatrix}).
}

\details{
  The \code{imdepi} data is a simplified version of what has been
  analyzed by Meyer et al. (2012). Simplification is with respect to the
  temporal resolution of the \code{stgrid} (see below) to be used in
  \code{\link{twinstim}}'s endemic model component.
  In what follows, we describe the elements \code{events},
  \code{stgrid}, \code{W}, and \code{qmatrix} of \code{imdepi} in
  greater detail.

  \code{imdepi$events} is a \code{"\linkS4class{SpatialPointsDataFrame}"}
  object (ETRS89 projection, i.e. EPSG code 3035, with unit \sQuote{km})
  containing 636 events, each with the following entries:
  \describe{
    \item{time:}{Time of the case occurrence measured in number of days
      since origin. Note that a U(0,1)-distributed random number has
      been subtracted from each of the original event times (days) to
      break ties (using
      \code{\link{untie}(imdepi_tied, amount=list(t=1))}).}
    \item{tile:}{Tile ID in the spatio-temporal grid (\code{stgrid}) of
      endemic covariates, where the event is contained in.
      This corresponds to one of the 413 districts of Germany.
    }
    \item{type:}{Event type, a factor with levels \code{"B"} and \code{"C"}.}
    \item{eps.t:}{Maximum temporal interaction range for the event.
      Here set to 30 days.}
    \item{eps.s:}{Maximum spatial interaction range for the event.
      Here set to 200 km.}
    \item{sex:}{Sex of the case, i.e. a factor with levels \code{"female"}
      and \code{"male"}. Note: for some cases this
      information is not available (\code{NA}).} 
    \item{agegrp:}{Factor giving the age group of the case,
      i.e. 0-2, 3-18 or >=19. Note: for one case this
      information is not available (\code{NA}).}
    \item{BLOCK, start:}{Block ID and start time (in days since origin) of
      the cell in the spatio-temporal endemic covariate grid, which the
      event belongs to.} 
    \item{popdensity:}{Population density (per square km) at the
      location of the event (corresponds to population density of the
      district where the event is located).}
  }
  There are further auxiliary columns attached to the events' data
  the names of which begin with a . (dot): These are created during
  conversion to the \code{"epidataCS"} class and are necessary for
  fitting the data with \code{twinstim}, see the description of the
  \code{"\link{epidataCS}"}-class.
  With \code{coordinates(imdepi$events)} one obtains the (x,y) locations
  of the events.

  The district identifier in \code{tile} is indexed according to
  the German official municipality key (
  \dQuote{Amtlicher Gemeindeschl\enc{ü}{ue}ssel}). See
  \url{https://de.wikipedia.org/wiki/Amtlicher_Gemeindeschl\%C3\%BCssel}
  for details.
  
  The data component \code{stgrid} contains the spatio-temporal grid of 
  endemic covariate information. In addition to the usual bookkeeping
  variables this includes:
  \describe{
    \item{area:}{Area of the district \code{tile} in square kilometers.}
    \item{popdensity:}{Population density (inhabitants per square
      kilometer) computed from DESTATIS (Federal Statistical Office)
      information (Date: 31.12.2008) on communities level (LAU2)
      aggregated to district level (NUTS3).} 
  }
  We have actually not included any time-dependent covariates here, we
  just established this grid with a (reduced -> fast) temporal
  resolution of \emph{monthly} intervals so that we can model endemic
  time trends and seasonality (in this discretized time).
  
  The entry \code{W} contains the observation window as a
  \code{"\linkS4class{SpatialPolygons}"} object, in this case the
  boundaries of Germany. It was obtained as
  \code{stateD <- rgeos::gUnaryUnion(districtsD)},
  where \code{districtsD} represents Germany's districts as at
  2009-01-01, simplified by the \dQuote{modified Visvalingam} algorithm
  (level 6.6\%) available at \url{https://MapShaper.org} (v. 0.1.17).
  The objects \code{districtsD} and \code{stateD} are contained in
  \code{system.file("shapes", "districtsD.RData", package="surveillance")}.
  
  The entry \code{qmatrix} is a \eqn{2\times 2}{2 x 2} identity matrix
  indicating that no transmission between the two finetypes can occur.
}

\source{
  IMD case reports: German Reference Centre for Meningococci (NRZM) --
  hosted by the Department of Hygiene and Microbiology,
  Julius-Maximilians-Universit\enc{ä}{ae}t W\enc{ü}{ue}rzburg, Germany.
  Thanks to Dr. Johannes Elias and Prof. Dr. Ulrich Vogel for providing
  the data. See \url{http://www.meningococcus.de/} and
  \url{http://www.episcangis.org/} for further details.

  Shapefile of Germany's districts as at 2009-01-01:
  German Federal Agency for Cartography and Geodesy, Frankfurt
  am Main, Germany, \url{https://gdz.bkg.bund.de/}.
  %% "Copy, distribution and making available to the public - also in
  %% parts - is allowed with reference."
}

\references{
  Meyer, S., Elias, J. and H\enc{ö}{oe}hle, M. (2012):
  A space-time conditional intensity model for invasive meningococcal
  disease occurrence. \emph{Biometrics}, \bold{68}, 607-616.
  \doi{10.1111/j.1541-0420.2011.01684.x}
}

\seealso{
  the data class \code{"\link{epidataCS}"}, and function
  \code{\link{twinstim}} for model fitting. 
}

\examples{
data("imdepi")

# Basic information
print(imdepi, n=5, digits=2)

# What is an epidataCS-object?
str(imdepi, max.level=4)
names(imdepi$events@data)
# => events data.frame has hidden columns
sapply(imdepi$events@data, class)
# marks and print methods ignore these auxiliary columns

# look at the B type only
imdepiB <- subset(imdepi, type == "B")
#<- subsetting applies to the 'events' component
imdepiB

# select only the last 10 events
tail(imdepi, n=10)   # there is also a corresponding 'head' method

# Access event marks
str(marks(imdepi))

# there is an update-method which assures that the object remains valid
# when changing parameters like eps.s, eps.t or qmatrix
update(imdepi, eps.t = 20)

# Summary
s <- summary(imdepi)
s
str(s)

# Step function of number of infectives
plot(s$counter, xlab = "Time [days]",
     ylab = "Number of infectious individuals",
     main = "Time series of IMD assuming 30 days infectious period")

# distribution of number of potential sources of infection
opar <- par(mfrow=c(1,2), las=1)
for (type in c("B","C")) {
  plot(100*prop.table(table(s$nSources[s$eventTypes==type])),
  xlim=range(s$nSources), xlab = "Number of potential epidemic sources",
  ylab = "Proportion of events [\%]")
}
par(opar)

# a histogram of the number of events along time (using the
# plot-method for the epidataCS-class, see ?plot.epidataCS)
opar <- par(mfrow = c(2,1))
plot(imdepi, "time", subset = type == "B", main = "Finetype B")
plot(imdepi, "time", subset = type == "C", main = "Finetype C")
par(opar)

# Plot the spatial distribution of the events in W
plot(imdepi, "space", points.args = list(col=c("indianred", "darkblue")),
     axes = TRUE, lwd = 2)
title(xlab = "x [km]", ylab = "y [km]")

\dontrun{
# or manually (no legends, no account for tied locations)
plot(imdepi$W, lwd=2)
plot(imdepi$events, pch=c(3,4)[imdepi$events$type], cex=0.8,
     col=c("indianred", "darkblue")[imdepi$events$type], add=TRUE)
}

\dontrun{
  # Show a dynamic illustration of the spatio-temporal dynamics of the 
  # spread during the first year of type B with a step size of 7 days
  animate(imdepiB, interval=c(0,365), time.spacing=7, sleep=0.1)
}
}

\keyword{datasets}
