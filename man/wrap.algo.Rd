\name{wrap.algo}
\alias{wrap.algo}
\alias{farrington}
\alias{bayes}
\alias{rki}
\alias{cusum}
\alias{glrpois}
\alias{glrnb}
\alias{outbreakP}
%% FIXME: hmm and rogerson are currently undocumented and unexported
%\alias{hmm}
%\alias{rogerson}
\encoding{latin1}

\title{Multivariate Surveillance through independent univariate algorithms}
\description{
  This function takes an \code{sts} object and applies an univariate
  surveillance algorithm to the time series of each observational unit.
}
\usage{
%This is the main function
wrap.algo(sts, algo, control,control.hook=function(k, control)
         return(control),verbose=TRUE,...)

%Derived functions fixing the control object and the "algo" argument
farrington(sts, control=list(range=NULL, b=5, w=3, reweight=TRUE,
         verbose=FALSE, alpha=0.05),...)
bayes(sts, control = list(range = range, b = 0, w = 6,
         actY = TRUE,alpha=0.05),...)
rki(sts, control = list(range = range, b = 2, w = 4,
         actY = FALSE),...)
cusum(sts,  control = list(range=range, k=1.04, h=2.26,
         m=NULL, trans="standard",alpha=NULL),...)
glrpois(sts, control = list(range=range,c.ARL=5, S=1,beta=NULL,
         Mtilde=1, M=-1, change="intercept",theta=NULL),...)
glrnb(sts, control = list(range=range,c.ARL=5, mu0=NULL, alpha=0,
         Mtilde=1, M=-1, change="intercept",
         theta=NULL,dir=c("inc","dec"),
         ret=c("cases","value")),...)
outbreakP(sts, control=list(range = range, k=100,
         ret=c("cases","value"),maxUpperboundCases=1e5),...)
}
\arguments{
  \item{sts}{Object of class \code{sts}}
  \item{algo}{Character string giving the function name of the algorithm
    to call, e.g. \code{"algo.farrington"}. Calling is done using
    \code{do.call}.}
  \item{control}{Control object as list. Depends on each algorithm.}
  \item{control.hook}{This is a function for handling multivariate
    objects. This argument is a function function of integer k and the
    current control object and which
    returns the appropriate control object for region k.}
  \item{verbose}{Boolean, if \code{TRUE} then textual information about the
    process is given}
  \item{...}{Additional arguments sent to the \code{algo} function.}
}
\value{
  An \code{sts} object with the \code{alarm}, \code{upperbound},
  etc. slots filled with the results of independent and univariate
  surveillance algorithm.
}
\seealso{
  \code{\link{algo.rki}}, \code{\link{algo.farrington}},
  \code{\link{algo.cusum}}, \code{\link{algo.glrpois}},
  \code{\link{algo.glrnb}}, \code{\link{algo.outbreakP}}
  for the exact form of the \code{control} object.
}
\author{M. \enc{Höhle}{Hoehle}}
\keyword{classif}
