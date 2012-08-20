\encoding{latin1}
\docType{data}
\name{imdepi}
\alias{imdepi}
\alias{imdepifit}

\title{
  Occurrence of Invasive Meningococcal Disease in Germany
}

\description{
  \code{imdepi} contains data on the spatio-temporal location of 636
  cases of invasive meningococcal disease caused by two specific
  meningococcal finetypes in Germany.
  \code{imdepifit} contains a model fit to the \code{imdepi} data. 
}

\usage{
data(imdepi)
data(imdepifit)
}

\format{
  \code{imdepi} is an object of class
  \code{"\link{epidataCS}"} (a list with components \code{events},
  \code{stgrid}, \code{W} and \code{qmatrix}).
  \code{imdepifit} is an object of class \code{"\link{twinstim}"}.
}

\details{
  The dataset contains a slightly modified version of the data analyzed in
  Meyer et al. (2012) represented as an \code{"epidataCS"} object, which
  is a list with components \code{events}, \code{stgrid}, \code{W}, and
  \code{qmatrix}.

  The \code{events} entry is a \code{"\linkS4class{SpatialPointsDataFrame}"}
  object (ETSR89 projection, i.e. EPSG code 3035, with unit \sQuote{km})
  containing 636 events, each with the following entries:
  \describe{
    \item{ID:}{Unique identifier of the case.}
    \item{time:}{Time of the case occurrence measured in number of days since origin.}
    \item{tile:}{Tile ID in the spatio-temporal grid (\code{stgrid}) of
      endemic covariates, where the event is contained in.
      This corresponds to one of the 413 districts of Germany.
    }
    \item{type:}{Event type, a factor with levels \code{"B"} and \code{"C"}.}
    \item{eps.t:}{Maximum temporal interaction range for the event.
      Here set to 30 days.}
    \item{eps.s:}{Maximum spatial interaction range for the event.
      Here set to 200 km.}
    \item{sex:}{Sex of the case, i.e. a factor with levels \code{"male"}
      and \code{"female"}.} 
    \item{agegrp:}{Factor giving the age group of the case,
      i.e. 0-2, 3-18 or >=19. Note: for a small number of cases this
      information is not available (\code{NA}).}
    \item{BLOCK, start:}{Block ID and start time (in days since origin) of
      the cell in the spatio-temporal endemic covariate grid, which the
      event belongs to.} 
    \item{popdensity:}{Population density at the location of the event
      (corresponds to population density of the district where the event is
      located).}
  }
  There are further auxiliary columns attached to the events' data
  the names of which begin with a . (dot): These are created during
  conversion to the \code{"epidataCS"} class and are necessary for
  fitting the data with \code{twinstim}, see the description of the
  \code{"\link{epidataCS}"}-class.
  With \code{coordinates(imdepi$events)} one obtains the (x,y) locations
  of the events.

  The identifier in \code{tile} is indexed according to
  the German official municipality key (
  \dQuote{Amtlicher Gemeindeschl\enc{ü}{ue}ssel}). See
  \url{http://de.wikipedia.org/wiki/Amtlicher_Gemeindeschl\%C3\%BCssel}
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
  resolution of monthly intervals so that we can model endemic time
  trends and seasonality (in this discretized time).
  
  The entry \code{W} contains the observation window as a
  \code{\link{SpatialPolygons}} object. In this case the boundaries of
  Germany (simplified polygon to enhance speed).
  %% The original shape file was downloaded from
  %% \url{http://www.geodatenzentrum.de/geodaten/gdz_rahmen.gdz_div?gdz_spr=eng&gdz_akt_zeile=5&gdz_anz_zeile=5&gdz_user_id=0}
  %% as of 2009-01-01, which has copyright by Bundesamt für Kartographie
  %% und Geodäsie, Frankfurt am Main:
  %% "Copy, distribution and making available to the public - also in
  %% parts - is allowed with reference."
  %% This shapefile was then simplified using \code{simplify.owin()} from
  %% package \pkg{spatstat} with parameter \code{dmin = 10}.
  
  The entry \code{qmatrix} is a \eqn{2\times 2}{2 x 2} identity matrix
  indicating that no transmission between the two finetypes can occur.

  The code used to produce the \code{imdepifit} can be found at the end
  of the examples section below.
}

\source{
  National Reference Centre for Meningococci (NRZM) -- hosted by the
  Department of Hygiene and Microbiology,
  Julius-Maximilians-Universit\enc{ä}{ae}t W\enc{ü}{ue}rzburg, Germany.
  Thanks to Dr. Johannes Elias and Prof. Dr. Ulrich Vogel from the NRZM
  for providing the IMD event data. See \url{http://www.meningococcus.de/}
  and \url{http://episcangis.hygiene.uni-wuerzburg.de/} for further
  details about the IMD event data.
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
plot(imdepi, aggregate = "time", subset = type == "B", main = "Finetype B")
plot(imdepi, aggregate = "time", subset = type == "C", main = "Finetype C")
par(opar)

# Plot the spatial distribution of the events in W
plot(imdepi$W)
colTypes <- c("darkblue","indianred2")
pch <- c(3,4)
with(imdepi, {
  plot(events,pch=pch[events$type],col=colTypes[events$type],add=TRUE)
  legend(x="topleft",levels(events$type),col=colTypes, pch=pch,title="Finetype")
})

# somehow easier: use the plot-method for the epidataCS class
plot(imdepi, aggregate="space", axes=TRUE, colTypes=colTypes)
title(xlab = "x [km]", ylab = "y [km]")

\dontrun{
  # Show a dynamic illustration of the spatio-temporal dynamics of the 
  # spread during the first year of type B with a step size of 7 days
  animate(imdepiB, interval=c(0,365), time.spacing=7, sleep=0.1)
}


\dontrun{
### data(imdepifit) is obtained as
data(imdepi)
imdepifit <- twinstim(
    endemic = ~1 + offset(log(popdensity)) + I(start/365) + 
        sin(start * 2 * pi/365) + cos(start * 2 * pi/365),
    epidemic = ~1 + type + agegrp, siaf = siaf.gaussian(1, logsd = TRUE),
    data = imdepi, subset = !is.na(agegrp), nCub = 24, nCub.adaptive = TRUE, 
    optim.args = list(par = c(-20, -0.041, 0.24, 0.35, -17, 0, 0, 0, 2.8)),
    model = FALSE, cumCIF = TRUE
)
}
}

\keyword{datasets}
