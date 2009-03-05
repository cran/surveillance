%This is just a way to "hide" functions, which the user should not know
%the proper way is to use NAMESPACE here (TODO!)
\name{surveillance-internal}
\alias{anscombeNB}
\alias{hypgeom2F1special}
\alias{hypGeomSmall}
%Some S4 stuff from sts
\alias{add.legend}
\alias{fix.dimnames}
\alias{hcl.colors}
\alias{insert.zeroes}
\alias{wait}
\alias{init.sts}
\alias{merge.list}
%Functions from M.Paul
\alias{commonDenom}
\alias{logInterpolation}
\alias{getH}
\alias{getHK}
\alias{getTheta1}
%additional (undocumented) functions needed for algo.hhh
\alias{unpackParams}
\alias{gradient}
\alias{weightedSumNeighbours}
\alias{jacobian}
\alias{sinCos2amplitudeShift}
\alias{jacobianAmplitudeShift}
\alias{getLambda}
\alias{expAlpha.mm}
\alias{logLik.ah}
\alias{logLik.ahg}
%
%twins function not documented
\alias{make.nu}
\alias{make.pois}
\alias{pois.plot}
\alias{tms.plot}

%Weird functions from sp package which produce warnings in R CMD check
%for the surveillance package (??)
\alias{summary,GridTopology-method}
\alias{summary,Spatial-method}

\title{Internal surveillance Functions}
\description{
  Undocumented internal \code{surveillance} functions. If you really
  want to know whats going on you might consider looking at the source
  code. If the package had a namespace you would not see this list. 
}
\details{
  These are not to be called by the user.
}
\keyword{ internal }
