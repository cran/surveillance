\name{stsNClist_animate}
\alias{stsNClist_animate}
\alias{animate_nowcasts}
\encoding{latin1}

\title{Animate a Sequence of Nowcasts}
\description{
  Animate a sequence of nowcasts stored as a list.
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
  \item{nowcasts}{A list of objects of class \code{\linkS4class{stsNC}}}
  \item{linelist_truth}{True linelist}
  \item{method}{Which method to show (has to be present in the
    nowcasts)}
  \item{control}{List with control options}
  \item{showLambda}{Boolean indicating whether to show the estimate for
    the epidemic curve (only applied to \code{bayes.trunc.ddcp})}
}
\value{
  This function is experimental and not yet fully documented.
}
\author{M. \enc{Höhle}{Hoehle}}
\seealso{
  \url{https://staff.math.su.se/hoehle/blog/2016/07/19/nowCast.html}
  for a worked through example.
}
\keyword{hplot}
