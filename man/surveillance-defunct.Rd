\name{surveillance-defunct}
\alias{surveillance-defunct}
\title{Defunct Functions in Package \pkg{surveillance}}
\alias{compMatrix.writeTable}
\alias{correct53to52}
\alias{enlargeData}
\alias{makePlot}
\alias{readData}
\alias{test}
\alias{testSim}
\alias{toFileDisProg}
\alias{algo.hhh}
\alias{algo.hhh.grid}
\alias{create.grid}
\description{
  The functions listed here are no longer part of \pkg{surveillance}.
}
\usage{
## Defunct in surveillance 1.18.0
algo.hhh(disProgObj, control=list(
             lambda=TRUE, neighbours=FALSE, linear=FALSE, nseason=0,
             negbin=c("none", "single", "multiple"),
             proportion=c("none", "single", "multiple"), lag.range=NULL
         ), thetastart=NULL, verbose=TRUE)
algo.hhh.grid(disProgObj, control=list(
                 lambda=TRUE, neighbours=FALSE, linear=FALSE, nseason=0,
                 negbin=c("none", "single", "multiple"),
                 proportion=c("none", "single", "multiple"), lag.range=NULL
              ), thetastartMatrix, maxTime=1800, verbose=FALSE)
create.grid(disProgObj, control, params=list(
                epidemic=c(0.1, 0.9, 5), endemic=c(-0.5, 0.5, 3),
                negbin=c(0.3, 12, 10)
            ))

## Defunct in surveillance 1.17.0
compMatrix.writeTable(compMatrix)
correct53to52(disProgObj, firstweek = 1)
enlargeData(disProgObj, range = 1:156, times = 1)
makePlot(outputpath, data = "k1", method = "rki1",
         name, disease, range = 157:339)
readData(abb, week53to52=TRUE, sysPath=TRUE)
test(data = c("k1", "m5"), range = 157:339)
testSim(p = 0.99, r = 0.01, length = 400, A = 1, alpha = 1,
        beta = 0, phi = 0, frequency = 1, state = NULL, K, 
        range = 200:400)
toFileDisProg(disProgObj, toFile)
}
\details{
  \code{algo.hhh} was an early implementation of the HHH regression
  framework for multivariate time series of infectious disease counts.
  An improved and considerably extended implementation is provided by
  the \code{\link{hhh4}} function since 2012.
  
  The trivial function \code{compMatrix.writeTable} is no longer used
  (it did nothing more than generating an \code{\link{xtable}} of its input).
  
  The ancient test function \code{makePlot} is no longer used,
  nor are functions \code{readData} (the datasets are still available
  via \code{data(\link{m1})} etc) and \code{correct53to52} /
  \code{enlargeData} (which both only worked for old \code{"disProg"}
  objects with non-matrix elements). \code{enlargeData} is still
  exemplified in the old \code{vignette("surveillance")}.
  \code{test} calls of outbreak detection algorithms applied to the old
  SurvStat datasets can be found in \code{vignette("surveillance")},
  and \code{testSim} is provided as an example
  in \code{help("\link{sim.pointSource}")}.

  Functions related to the old \code{"\link[=create.disProg]{disProg}"}
  class are no longer needed. The effect of \code{toFileDisProg} could
  still be achieved via \code{\link{write.table}} of
  \code{as.data.frame(disProg2sts(disProgObj))[c("epoch", "observed", "state")]}.
}
\seealso{
  \code{\link{Defunct}}
}
\keyword{internal}
