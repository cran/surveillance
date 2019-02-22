\encoding{latin1}
\name{hagelloch}
\alias{hagelloch}
\alias{hagelloch.df}
\docType{data}
\keyword{datasets}

\title{1861 Measles Epidemic in the City of Hagelloch, Germany}
\description{
  Data on the 188 cases in the measles outbreak among children in the
  German city of Hagelloch (near T\enc{ü}{ue}bingen) 1861. The data were
  originally collected by Dr. Albert Pfeilsticker (1863) and augmented and
  re-analysed by Dr. Heike Oesterle (1992).
  This dataset is used to illustrate the \code{twinSIR} model class in
  \code{vignette("twinSIR")}.
}

\usage{
data("hagelloch")
}

\format{
  Loading \code{data("hagelloch")} gives two objects:
  \code{hagelloch} and \code{hagelloch.df}.
  The latter is the original \code{data.frame} of 188 rows
  with individual information for each infected child.
  \code{hagelloch} has been generated from \code{hagelloch.df}
  via \code{\link{as.epidata}} (see the Examples below) to obtain an
  \code{"epidata"} object for use with \code{\link{twinSIR}}.
  It contains the entire SIR event history of the outbreak
  (but not all of the covariates).

  The covariate information in \code{hagelloch.df} is as follows:
  
  \describe{
    \item{PN:}{patient number}
    \item{NAME:}{patient name (as a factor)}
    \item{FN:}{family index}
    \item{HN:}{house number}
    \item{AGE:}{age in years}
    \item{SEX:}{gender of the individual (factor: male, female)}
    \item{PRO:}{\code{Date} of prodromes}
    \item{ERU:}{\code{Date} of rash}
    \item{CL:}{class (factor: preschool, 1st class, 2nd class)}
    \item{DEAD:}{\code{Date} of death (with missings)}
    \item{IFTO:}{number of patient who is the putative source of infection (0 = unknown)}
    \item{SI:}{serial interval = number of days between dates of prodromes of infection source and infected person}
    \item{C:}{complications (factor: no complications, bronchopneumonia, severe bronchitis, lobar pneumonia, pseudocroup, cerebral edema)}
    \item{PR:}{duration of prodromes in days}
    \item{CA:}{number of cases in family}
    \item{NI:}{number of initial cases}
    \item{GE:}{generation number of the case}
    \item{TD:}{day of max. fever (days after rush)}
    \item{TM:}{max. fever (degree Celsius)}
    \item{x.loc:}{x coordinate of house (in meters). Scaling in metres
      is obtained by multiplying the original coordinates by 2.5 (see details
      in Neal and Roberts (2004))}
    \item{y.loc:}{y coordinate of house (in meters). See also the above
      description of \code{x.loc}.}
    \item{tPRO:}{Time of prodromes (first symptoms) in days after the start of the
      epidemic (30 Oct 1861).}
    \item{tERU:}{Time upon which the rash first appears.}
    \item{tDEAD:}{Time of death, if available.}
    \item{tR:}{Time at which the infectious period of the individual is
      assumed to end. This unknown time is calculated as
      \deqn{tR_i = \min{tDEAD_i,tERU_i+d_0},}{tR[i] = min(tDEAD[i],tERU[i]+d0),}
      where -- as in Section 3.1 of Neal and Roberts (2004) -- we use
      \eqn{d_0=3}{d0=3}.}
    \item{tI:}{Time at which the individual is assumed to become
      infectious. Actually this time is unknown, but we use
      \deqn{tI_i = tPRO_i - d_1,}{tI[i] = tPRO[i] - d1,} 
      where \eqn{d_1=1}{d1=1} as in Neal and Roberts (2004).}
  }

  The time variables describe the transitions of the individual in an
  Susceptible-Infectious-Recovered (SIR) model.  Note that in
  order to avoid ties in the event times resulting from daily
  interval censoring, the times have been jittered uniformly within the
  respective day. The time point 0.5 would correspond to noon of 30 Oct 1861.

  The \code{hagelloch} \code{"epidata"} object only retains some of
  the above covariates to save space. Apart from the usual
  \code{"epidata"} event columns, \code{hagelloch} contains a number of
  extra variables representing distance- and covariate-based weights for
  the force of infection:
  \describe{
    \item{household:}{the number of currently infectious children in the same
      household (including the child itself if it is currently infectious).}
    \item{nothousehold:}{the number of currently infectious children
      outside the household.}
    \item{c1, c2:}{the number of children infectious during the respective
      time block and being members of class 1 and 2, respectively; but
      the value is 0 if the individual of the row is not herself a
      member of the respective class.}
  }
  Such epidemic covariates can been computed by specifying suitable
  \code{f} and \code{w} arguments in \code{\link{as.epidata}} at
  conversion (see the code below), or at a later step via the
  \code{\link[=update.epidata]{update}}-method for \code{"epidata"}.
}

\source{
  Thanks to Peter J. Neal, University of Manchester, for providing us
  with these data, which he again became from Niels Becker, Australian
  National University. To cite the data, the main references are Pfeilsticker
  (1863) and Oesterle (1992). 
}

\examples{
data("hagelloch")
head(hagelloch.df)   # original data documented in Oesterle (1992)
head(as.data.frame(hagelloch))   # "epidata" event history format

## How the "epidata" 'hagelloch' was created from 'hagelloch.df'
stopifnot(all.equal(hagelloch,
  as.epidata(
    hagelloch.df, t0 = 0, tI.col = "tI", tR.col = "tR",
    id.col = "PN", coords.cols = c("x.loc", "y.loc"),
    f = list(
        household    = function(u) u == 0,
        nothousehold = function(u) u > 0
    ),
    w = list(
        c1 = function (CL.i, CL.j) CL.i == "1st class" & CL.j == CL.i,
        c2 = function (CL.i, CL.j) CL.i == "2nd class" & CL.j == CL.i
    ),
    keep.cols = c("SEX", "AGE", "CL"))
))


### Basic plots produced from hagelloch.df

# Show case locations as in Neal & Roberts (different scaling) using
# the data.frame (promoted to a SpatialPointsDataFrame)
coordinates(hagelloch.df) <- c("x.loc","y.loc")
plot(hagelloch.df, xlab="x [m]", ylab="x [m]", pch=15, axes=TRUE,
     cex=sqrt(multiplicity(hagelloch.df)))

# Epicurve
hist(as.numeric(hagelloch.df$tI), xlab="Time (days)", ylab="Cases", main="")


### "epidata" summary and plot methods

(s <- summary(hagelloch))
head(s$byID)
plot(s)

\dontrun{
  # Show a dynamic illustration of the spread of the infection
  animate(hagelloch, time.spacing=0.1, sleep=1/100,
          legend.opts=list(x="topleft"))
}
}

\references{
  Pfeilsticker, A. (1863). Beitr\enc{ä}{ae}ge zur Pathologie der Masern mit
  besonderer Ber\enc{ü}{ue}cksichtigung der statistischen Verh\enc{ä}{ae}ltnisse,
  M.D. Thesis, Eberhard-Karls-Universit\enc{ä}{ae}t T\enc{ü}{ue}bingen.
  Available as \url{http://www.archive.org/details/beitrgezurpatho00pfeigoog}.
  
  Oesterle, H. (1992). Statistische Reanalyse einer Masernepidemie 1861
  in Hagelloch, M.D. Thesis, Eberhard-Karls-Universit\enc{ä}{ae}at
  T\enc{ü}{ue}bingen.

  Neal, P. J. and Roberts, G. O (2004). Statistical inference and model
  selection for the 1861 Hagelloch measles epidemic, Biostatistics
  5(2):249-261
}

\seealso{
  data class: \code{\link{epidata}}
  
  point process model: \code{\link{twinSIR}}
  
  illustration with \code{hagelloch}: \code{vignette("twinSIR")}
}
