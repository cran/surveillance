\name{stsNClist_animate}
\alias{stsNClist_animate}
\alias{animate_nowcasts}
\encoding{latin1}

\title{Animate a sequence of nowcasts}
\description{Animate a sequence of nowcasts stored as a list.
}
\usage{
animate_nowcasts(nowcasts,linelist_truth, method="bayes.trunc.ddcp",
      control=list(dRange=NULL,anim.dRange=NULL, plot.dRange=NULL,
                   consistent=FALSE, sys.sleep = 1, ylim=NULL,cex.names=0.7,
                   col=c("violetred3","#2171B5","orange","blue","black",
                         "greenyellow")),
                   showLambda=TRUE)
}
\arguments{
  \item{nowcasts}{A list of objects of class \code{stsNC}}
  \item{linelist_truth}{True linelist}
  \item{method}{Which method to show (has to be present in the
    nowcasts)}
  \item{control}{List with control options}
  \item{showLambda}{Boolean indicating whether to show the estimate for
    the epidemic curve (only applied to \code{bayes.trunc.ddcp})}
}
\value{
  This function is experimental and is not yet documented.
}
\details{
  This function is experimental and might be changed in the future.
}
\author{M. \enc{Höhle}{Hoehle}}
\examples{
## See http://staff.math.su.se/hoehle/blog/2016/07/19/nowCast.html for
## a worked through example. Code will migrate into the package in due
## course.
}
\keyword{hplot}

