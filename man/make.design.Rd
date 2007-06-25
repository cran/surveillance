\name{make.design}
\alias{make.design}
\title{Create the design matrices}
\description{
        Creates the design matrices needed for \code{meanResponse}
}
\usage{
  make.design(disProgObj, control=list(linear=FALSE, nseason=0, period=52,
      neighbours=FALSE, negbin=FALSE, lambda=TRUE))
}
\arguments{
        \item{disProgObj}{Object of class \code{disProg}}
        \item{control}{Control object:
              \itemize{
              \item{linear}{if \code{true} a linear trend is included.}
              \item{nseason}{number of Fourier frequencies to be used}
              \item{period}{period length (defaults to 52 for weekly data)}
              \item{neighbours}{if \code{true} a autoregressive parameter for 
                    adjacent areas is included.}
              \item{negbin}{if \code{true} negative binomial rather than poisson is used.}
              \item{lambda}{if \code{true} a autoregressive parameter is included.}
              }}
}
\value{\item{list}{
    \itemize{
    \item{Y}{matrix with number of cases \eqn{y_{it}} in unit \eqn{i} at 
          time \eqn{t} as elements, i.e. data without the first time point.} 
    \item{Ym1}{matrix with previous number of cases \eqn{y_{i,t-1}},
          i.e data without the last time point.}
    \item{Ym1.neighbours}{matrix with sum of earlier counts of adjacent units
          \eqn{\sum_{j \sim i} y_{j,t-1}} }
    \item{X.trendSeason}{design matrix for linear trend and seasonal components}
    \item{X.predict}{design matrix for prediction  } %
    \item{populationFrac}{matrix with corresponding population proportions} 
    \item{dimtheta}{number of parameters used in model}
    \item{psiIndex}{position of \eqn{\psi} in vector with parameters}
    \item{control}{control object}
    }}
}
\details{}
\seealso{}
\author{M.Paul, L. Held}
\examples{}
\keyword{}
\source{}

