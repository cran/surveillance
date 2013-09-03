\name{stcd}
\alias{stcd}
\encoding{latin1}

\title{Spatio-temporal cluster detection}
\description{
  Shiryaev-Roberts based prospective spatio-temporal cluster detection
  as in Assuncao & Correa (2009).
}
\usage{
stcd(x, y,t,radius,epsilon,areaA, areaAcapBk, threshold, cusum=FALSE)
}
\arguments{
  \item{x}{Vector containing spatial x coordinate of the events.}
  \item{y}{Vector containing spatial y coordinate of the events.}
  \item{t}{Vector containing the time points of the events. It is
    assumed that the vector is sorted (early->last).}
  \item{radius}{Radius of the cluster to detect.}
  \item{epsilon}{Relative change of event-intensity within the cluster
    to detect. See reference paper for an explicit definition.}
  \item{areaA}{Area of the observation region A (single number) -- This
    argument is currently ignored!}
  \item{areaAcapBk}{Area of A \ B(s_k,rho) for all k=1,\ldots,n
    (vector). This argument is currently ignored!}
  \item{threshold}{Threshold limit for the alarm and should be equal to
    the desired Average-Run-Length (ARL) of the detector.}
  \item{cusum}{(logical) If \code{FALSE} (default) then the
    Shiryaev-Roberts detector is used as in the original article by
    Assuncao & Correa (2009), i.e. \eqn{R_n = \sum_{k=1}^n
      \Lambda_{k,n}}, where \eqn{\Lambda_{k,n}} denotes the likelihood
    ratio between the in-control and out-of control model. If
    \code{TRUE}, CUSUM test statistic is 
    used instead. Here, \deqn{R_n = \max_{1\leq k \leq n}
      \Lambda_{k,n}}. Note that this has implications on what threshold
    will sound the alarm (CUSUM threshold needs to be smaller).}
}
\details{
 Shiryaev-Roberts based spatio-temporal cluster detection based
 on the work in Assuncao and Correa (2009). The implementation
 is based on C++ code originally written by Marcos Oliveira Prates, UFMG,
 Brazil and provided by Thais Correa, UFMG, Brazil during her research
 stay in Munich. This stay was financially supported by the Munich
 Center of Health Sciences.

 Note that the vectors \code{x}, \code{y} and \code{t} need to be of the
 same length. Furthermore, the vector \code{t} needs to be sorted (to
 improve speed, the latter is not verified within the function).
 
 The current implementation uses a call to a C++ function to perform the
 actual computations of the test statistic. The function is currently
 experimental -- data type and results may be subject to changes.

}
\value{A list with three components
  \item{R}{A vector of the same length as the input containing the value
    of the test statistic for each observation.}
  \item{idxFA}{Index in the x,y,t vector causing a possible alarm. If no
    cluster was detected, then a value of \code{-1} is returned here.}
  \item{idxCC}{index in the x,y,t vector of the event containing the
    cluster. If no cluster was detected, then a value of \code{-1} is
    returned here.}
}
\references{
Assuncao, R. and Correa, T. (2009), Surveillance to detect emerging
space-time clusters, Computational Statistics & Data Analysis,
53(8):2817-2830.
}
\examples{
if (require("splancs")) {
    # load the data from package "splancs"
    data(burkitt, package="splancs")

    # order the times
    burkitt <- burkitt[order(burkitt$t), ]

    #Parameters for the SR detection
    epsilon <- 0.5 # relative change within the cluster
    radius <- 20 # radius
    threshold <- 161 # threshold limit

    res <- stcd(x=burkitt$x,
                y=burkitt$y,
                t=burkitt$t,
                radius=radius,
                epsilon=epsilon,
                areaA=1,
                areaAcapBk=1,
                threshold=threshold)

    #Index of the event
    which.max(res$R >= threshold)
}
}
\author{M. O. Prates, T. Correa and M. \enc{Höhle}{Hoehle}}
\keyword{cluster}

