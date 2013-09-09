\name{hagelloch}
\alias{hagelloch}
\alias{hagelloch.df}
\encoding{latin1}

\title{1861 measles epidemic in the city of Hagelloch, Germany}
\description{Data on the 188 cases in the measles outbreak among children in the
  German city of Hagelloch (near T\enc{ü}{ue}bingen) 1861. The data were
  originally collected by Dr. Albert Pfeilsticker and re-analysed by
  Dr. Heike Oesterle.  
}
\usage{
data("hagelloch")
}
\format{
  Loading the data in \code{hagelloch} gives two objects:
  \code{hagelloch} and \code{hagelloch.df}. The former is an
  \code{epidata} object for use with \code{twinSIR} containing the
  entire process history of outbreak in the population of 188 children.
  The latter is a simple \code{data.frame} containing a single line
  with individual information of each infected child. 

  The covariate information in either \code{hagelloch} or
  \code{hagelloch.df} are as follows:
  
  \describe{
    \item{PN}{patient no.}
    \item{NAME}{patient name}
    \item{FN}{family name}
    \item{HN}{House number.}
    \item{AGE}{Age in years}
    \item{SEX}{Sex of the individual (unknown,male,female)}
    \item{PRO}{Date of prodromes, formatted as \code{Date}}
    \item{ERU}{Date of rash, formatted as \code{Date}}
    \item{CL}{Class (preschool, 1st class, 2nd class)}
    \item{DEAD}{Date of death (as \code{Date}), if not available then \code{NA}}
    \item{IFTO}{no. of patient who most likely insource of infection (0 = unknown)}
    \item{SI}{serialinterval = no of days between dates of prodromes between infection source and infected person}
    \item{C}{complications (no complicatons, bronchopneumonia, severe bronchitis, lobar pneumonia, pseudocroup, cerebral edema)}
    \item{PR}{duration of prodromes in days}
    \item{CA}{number of cases in family}
    \item{NI}{number of initial cases}
    \item{GE}{generation number of case}
    \item{TD}{day of max. fever (days after rush)}
    \item{TM}{man fever (grad celsius)}
    \item{HNX}{x coordinate of house (in meters). Scaling in metres
      is obtained by taking original data multiplied by 2.5 (see details
      in Neal and Roberts (2004))}
    \item{HNY}{y coordinate of house (in meters). See also the above
      description of \code{HNX}.}
    \item{tS}{Time of symptoms (here: prodomes) in days after the start of the
      epidemic 30 Oct 1861. This corresponds to the time of the
    "Exposed"->"Infectious" transition in an SEIR model.}
    \item{tQ}{Time upon which the rash first appears.}
    \item{tD}{Time of death, if available, in day after start of the
      epidemic.} 
    \item{tI}{Time at which individual becomes infected. This
      corresponds to the transition from state "Susceptible" to
      "Exposed" in an "SEIR" model. Actually
      this time is unknown, but we use
      \deqn{I_i = S_i - d_1,}{I[i] = S[i] - d[1],} where  
      \eqn{d_1=1}{d[1]=1} is the latency period. See Neal and Roberts
    (2004) for details.} 
    \item{tR}{Time at which the infectious period of the individual ends,
      i.e. time at which the individual moves from "Infectious" to
      "Recovered". This unknown
      time is calculated as
      \deqn{R_i = \min{D_i,Q_i+d_0},}{R[i] = min(D[i],Q[i]+d[0]),}
      where as in Section 3.1 of Neal and Roberts (2004) we use
      \eqn{d_0=3}{d0=3}.} 
  }

  Note that the \code{tS},\code{tQ},\code{tD},\code{tI} and \code{tR}
  entries are only available in the \code{hagelloch} object. They
  describe the transition of individual in an
  Susceptible-Exposed-Infectious-Revovered (SEIR) model.  Note that in
  order to avoid ties in the event times resulting from the daily
  interval censoring, a \eqn{U(0,1)} distributed random variable has
  been subtracted from \code{tS}, \code{tD} and \code{tQ} in order to
  break ties.
  
}

\source{
  Thanks to Peter J. Neal, University of Manchester, for providing us
  with these data, which he again became from Niels Becker, Australian
  National University. To cite the data, the main references are Pfeilsticker
  (1863) and Oesterle (1992). 
}

\examples{
  data("hagelloch")
  #Data frame with the original data documented in Oesterle (1992)
  head(hagelloch.df)
  #Represented as epidata object
  (s <- summary(hagelloch))

  ### Use data.frame to illustrate the basics

  #Show case locations as in Neal & Roberts (different scaling) using
  #the data.frame (promoted to a SpatialPointsDataFrame)
  coordinates(hagelloch.df) <- c("HNX","HNY")
  plot(hagelloch.df,xlab="x [m]",ylab="x [m]",pch=15,axes=TRUE,
       cex=sqrt(multiplicity(hagelloch.df)))

  #Epicurve
  hist(as.numeric(hagelloch.df$tS),xlab="Time (days)",ylab="Cases",main="")

  ### SIR model information for population & individuals

  plot(s, col=c("green","red","darkgray"))
  stateplot(s, id=c("187"))

  \dontrun{
    #Show a dynamic illustration of the spread of the infection
    animate(hagelloch,time.spacing=0.1,legend.opts=list(x="topleft"),sleep=1/100)
  }
}

\references{
  Pfeilsticker, A. (1863). Beitr\enc{ä}{ae}ge zur Pathologie der Masern mit
  besonderer Ber\enc{ü}{ue}cksichtgung der statistischen Verh\enc{ä}{ae}ltnisse,
  M.D. Thesis, Eberhard-Karls Universit\enc{ä}{ae}t, T\enc{ü}{ue}bingen.
  Available as \url{http://www.archive.org/details/beitrgezurpatho00pfeigoog}.
  
  Oesterle, H. (1992). Statistische Reanalyse einer Masernepidemie 1861
  in Hagelloch, M.D. Thesis, Eberhard-Karls Universit\enc{ä}{ae}at,
  T\enc{ü}{ue}bingen.

  Neal, P. J. and Roberts, G. O (2004). Statistical inference and model
  selection for the 1861 Hagelloch measles epidemic, Biostatistics
  5(2):249-261
}
\seealso{\code{\link{twinSIR}}, \code{\link{epidata}}}
\keyword{datasets}

