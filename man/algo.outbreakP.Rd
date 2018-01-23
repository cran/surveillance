\encoding{latin1}
\name{algo.outbreakP}
\alias{algo.outbreakP}
\alias{calc.outbreakP.statistic}
\title{Semiparametric surveillance of outbreaks}
\description{
  Frisen and Andersson (2009) method for semiparametric surveillance of outbreaks
}
\usage{
algo.outbreakP(disProgObj, control = list(range = range, k=100,
               ret=c("cases","value"),maxUpperboundCases=1e5))
}
\arguments{
\item{disProgObj}{object of class disProg (including the observed and the state chain).}
\item{control}{A list controlling the behaviour of the algorithm
    \describe{
    \item{\code{range}}{determines the desired
        time-points which should be monitored. Note that it is
        automatically assumed that ALL other values in \code{disProgObj}
        can be used for the estimation, i.e. for a specific value \code{i}
        in \code{range} all values from 1 to \code{i} are used for estimation.}
    \item{\code{k}}{The threshold value. Once the outbreak statistic
        is above this threshold \code{k} an alarm is sounded.}
    \item{\code{ret}}{a string specifying the type of
        \code{upperbound}-statistic that is returned. With
        \code{"cases"} the number of cases that would have been
        necessary to produce an alarm (NNBA) or with \code{"value"} the
        outbreakP-statistic is computed (see below).}
    \item{\code{maxUpperboundCases}}{Upperbound when numerically searching
        for NNBA. Default is 1e5.}
    }
  }
}
\value{
  \code{algo.outbreakP} gives a list of class \code{survRes} which
  includes the vector of alarm values for every time-point in
  \code{range}, the vector of threshold values for every time-point
  in \code{range}.
}
\details{
  A generalized likelihood ratio test based on the Poisson
  distribution is implemented where the means of the in-control and
  out-of-control states are computed by isotonic regression.

 \deqn{OutbreakP(s) = \prod_{t=1}^s \left( \frac{\hat{\mu}^{C1}(t)}{\hat{\mu}^D(t)} \right)^{x(t)}}

 where \eqn{\hat{\mu}^{C1}(t)} is the estimated mean obtained by
 uni-modal regression under the assumption of one change-point and
 \eqn{\hat{\mu}^D(t)} is the estimated result when there is no
 change-point (i.e. this is just the mean of all observations). Note
 that the contrasted hypothesis assume all means are equal until the
 change-point, i.e. this detection method is especially suited for
 detecting a shift from a relative constant mean. Hence, this is less
 suited for detection in diseases with strong seasonal endemic
 component. Onset of influenza detection is an example where this
 method works particular well.

 In case \code{control$ret == "cases"} then a brute force numerical
 search for the number needed before alarm (NNBA) is performed. That
 is, given the past observations, whats the minimum number which would
 have caused an alarm? Note: Computing this might take a while because
 the search is done by sequentially increasing/decreasing the last
 observation by one for each time point in \code{control$range} and
 then calling the workhorse function of the algorithm again. The argument
 \code{control$maxUpperboundCases} controls the upper limit of this
 search (default is 1e5).
 Currently, even though the statistic has passed the threshold, the NNBA
 is still computed. After a few time instances what typically happens is
 that no matter the observed value we would have an alarm at this time point. In this case the value of NNBA is set to \code{NA}. Furthermore, the first time
 point is always \code{NA}, unless \code{k<1}.

}
\source{
 The code is an extended R port of the Java code by Marianne
 \enc{Frisén}{Frisen} and Linus \enc{Schiöler}{Schioeler} from the
 CASE project available under the GNU GPL License v3. See
 \url{https://case.folkhalsomyndigheten.se/} for further details on the CASE
 project. A manual on how to use an Excel implementation of the method
 is available at \url{http://economics.handels.gu.se/english/Units+and+Centra/statistical_research_unit/software}.

 An additional feature of the R code is that it contains a search for
 NNBA (see details).
}
\author{M. \enc{Höhle}{Hoehle} -- based on Java code by M. Frisen and
  L. \enc{Schiöler}{Schioeler}}
\references{
\enc{Frisén}{Frisen}, M., Andersson and \enc{Schiöler}{Schioeler}, L., (2009), Robust
outbreak surveillance of epidemics in Sweden, Statistics in
Medicine, 28(3):476-493.

\enc{Frisén}{Frisen}, M. and Andersson, E., (2009) Semiparametric
Surveillance of Monotonic Changes, Sequential Analysis 28(4):434-454.
}
\examples{
#Use data from outbreakP manual (http://www.hgu.gu.se/item.aspx?id=16857)
y <- matrix(c(1,0,3,1,2,3,5,4,7,3,5,8,16,23,33,34,48),ncol=1)

#Generate sts object with these observations
mysts <- sts(y, alarm=y*0)

#Run the algorithm and present results
#Only the value of outbreakP statistic
upperbound(outbreakP(mysts, control=list(range=1:length(y),k=100,
           ret="value")))

#Graphical illustration with number-needed-before-alarm (NNBA) upperbound.
res <- outbreakP(mysts, control=list(range=1:length(y),k=100,
           ret="cases"))
plot(res,dx.upperbound=0,lwd=c(1,1,3),legend.opts=list(legend=c("Infected",
      "NNBA","Outbreak","Alarm"),horiz=TRUE))
}
\keyword{classif}
