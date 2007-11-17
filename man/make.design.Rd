\name{make.design}
\alias{make.design}
\encoding{latin1}
\title{Create the design matrices}
\description{
        Creates the design matrices needed for \code{meanResponse}
}
\usage{
make.design(disProgObj, control=list(lambda=TRUE, neighbours=FALSE, 
   linear=FALSE, nseason=0,
   negbin=c("none", "single", "multiple"), 
   proportion=c("none", "single", "multiple")) )
}
\arguments{
        \item{disProgObj}{Object of class \code{disProg}}
        \item{control}{Control object:
              \itemize{
              \item{lambda}{if \code{TRUE} an autoregressive 
                   parameter \eqn{\lambda} is included, if \code{lambda} is a vector of logicals,
                   unit-specific parameters \eqn{\lambda_i} are included}
              \item{neighbours}{if \code{TRUE} an autoregressive parameter for 
                    adjacent units \eqn{\phi} is included, if \code{neighbours} is a vector of logicals,
                    unit-specific parameters \eqn{\phi_i} are included}
              \item{linear}{a \code{logical} (or a vector of logicals) indicating wether a linear 
                    trend \eqn{\beta} (or a linear trend \eqn{\beta_i} for each unit)
                    is included}
              \item{nseason}{integer number of Fourier frequencies \code{nseason}; if \code{nseason} is a vector
                    of integers, each unit \eqn{i} gets its own seasonal parameters }
              \item{negbin}{if \code{"single"} negative binomial rather than poisson is used,
                    if \code{"multiple"} unit-specific overdispersion parameters are used.}
              \item{proportion}{see details in \code{\link{meanResponse}} }
              }}
}
\value{\item{list}{
    \itemize{
    \item{Y}{matrix with number of cases \eqn{y_{it}}{y_it} in unit \eqn{i} at 
          time \eqn{t} as elements, i.e. data without the first time point.} 
    \item{Ym1}{matrix with previous number of cases \eqn{y_{i,t-1}}{y_i,t-1},
          i.e data without the last time point.}
    \item{Ym1.neighbours}{matrix with weighted sum of earlier counts of adjacent units
          \eqn{\sum_{j \sim i} m_{ji} y_{j,t-1}}{sum_(j ~ i) w_ji * y_j,t-1} }
    \item{nOfNeighbours}{vector with number of neighbours for each unit \eqn{i} }
    \item{X.trendSeason}{design matrix for linear trend and seasonal components}
    \item{populationFrac}{matrix with corresponding population proportions} 
    \item{dimTheta}{list with number of parameters used in model}
    \item{control}{control object}
    \item{disProgObj}{Object of class \code{disProg}}
    }}
}
\author{M.Paul, L. Held}
\keyword{misc}



% predict.ah
