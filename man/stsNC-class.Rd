\name{stsNC-class}
\docType{class}
\alias{stsNC-class}
%New stsNC specific methods
\alias{reportingTriangle}
\alias{reportingTriangle,stsNC-method}
\alias{delayCDF}
\alias{delayCDF,stsNC-method}
\alias{score}
\alias{score,stsNC-method}
\alias{predint}
\alias{predint,stsNC-method}
%Coerce method to convert to sts object
\alias{coerce,sts,stsNC-method}

\encoding{latin1}

\title{Class "stsNC" -- a class inheriting from class \code{sts} which
  allows the user to store the results of back-projecting
  surveillance time series}

\description{
  A class inheriting from class \code{sts}, but with additional slots
  to store the results of nowcasting.
}

\section{Slots}{
  The slots are as for \code{"\linkS4class{sts}"}. However, a number of
  additional slots exists.

  \describe{
    \item{\code{reportingTriangle}:}{An array containing the upper and lower limit
      of the confidence interval.}
    \item{\code{predPMF}:}{Predictive distribution for each nowcasted
      time point.}
    \item{\code{pi}:}{A prediction interval for each nowcasted time
      point. This is calculated based on \code{predPMF}.}
    \item{\code{truth}:}{An object of type \code{sts} containing the
      true number of cases.}
    \item{\code{delayCDF}:}{List with the CDF of the estimated delay
      distribution for each method.}
    \item{\code{SR}:}{Possible output of proper scoring rules}
  }
}
\section{Methods}{
  The methods are the same as for \code{"\linkS4class{sts}"}.

  \itemize{
    \item{\code{signature(from = "sts", to = "stsNC")}}{
      Convert an object of class \code{sts} to class \code{stsNC}.
    }
    \item{reportingTriangle}{\code{signature(x = "stsNC")}: extract the
      \code{reportingTriangle} slot of an \code{stsNC} object.
    }
    \item{delayCDF}{\code{signature(x = "stsNC")}: extract the
      \code{delayCDF} slot of an \code{stsNC} object.
    }
    \item{score}{\code{signature(x = "stsNC")}: extract the
      scoring rules result slot of an \code{stsNC} object.
    }
    \item{predint}{\code{signature(x = "stsNC")}: extract the
      prediction interval slot of an \code{stsNC} object.
    }
}
}

\author{M. \enc{Höhle}{Hoehle}}

\keyword{classes}
