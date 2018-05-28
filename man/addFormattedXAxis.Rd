\encoding{latin1}
\name{addFormattedXAxis}
\alias{addFormattedXAxis}
% helper functions for time axis formatting
\alias{atChange}
\alias{at2ndChange}
\alias{atMedian}

\title{
Formatted Time Axis for \code{"sts"} Objects
}

\description{
Add a nicely formatted x-axis to time series plots related to the
\code{"\linkS4class{sts}"} class. This utility function is, e.g., used
by \code{\link{stsplot_time1}} and \code{\link{plotHHH4_fitted1}}.
}

\usage{
addFormattedXAxis(x, epochsAsDate = FALSE,
                  xaxis.tickFreq = list("\%Q"=atChange),
                  xaxis.labelFreq = xaxis.tickFreq,
                  xaxis.labelFormat = "\%G\n\n\%OQ",
                  ...)
}

\arguments{
  \item{x}{
    an object of class \code{"\linkS4class{sts}"}.
  }
  \item{epochsAsDate}{
    a logical indicating if the old (\code{FALSE}) or the new
    (\code{TRUE}) and more flexible implementation should be used.
    The \code{xaxis.*} arguments are only relevant for the new
    implementation \code{epochsAsDate = TRUE}.
  }
  \item{xaxis.labelFormat,xaxis.tickFreq,xaxis.labelFreq}{
    see the details below.
  }
  \item{\dots}{
    further arguments passed to \code{\link{axis}}.
  }
}

\details{
  The setting \code{epochsAsDate = TRUE}
  enables very flexible formatting of the x-axis and its
  annotations using the \code{xaxis.tickFreq}, \code{xaxis.labelFreq}
  and \code{xaxis.labelFormat} arguments. The first two are named lists containing
  pairs with the \emph{name} being a \code{\link{strftime}} single
  conversion specification and the second part is a function which based
  on this conversion returns a subset of the rows in the \code{sts}
  objects. The subsetting function has the following header:
  \code{function(x,xm1)}, where \code{x} is a vector containing
  the result of applying the conversion in \code{name} to the epochs of
  the \code{sts} object and \code{xm1} is the scalar result when
  applying the conversion to the natural element just before the first
  epoch. Please note that the input to the subsetting function is converted
  using \code{as.numeric} before calling the function. Hence, the
  conversion specification needs to result in a string convertible to integer.
  
  Three predefined subsetting functions exist:
  \code{atChange}, \code{at2ndChange}  and \code{atMedian}, which
  are used to make a tick at each (each 2nd for \code{at2ndChange})
  change and at the median index computed on all having the same value,
  respectively:
  \preformatted{
    atChange <- function(x,xm1) which(diff(c(xm1,x)) != 0)
    at2ndChange <- function(x,xm1) which(diff(c(xm1,x) \%/\% 2) != 0)
    atMedian <- function(x,xm1) tapply(seq_along(x), INDEX=x, quantile, prob=0.5, type=3)
  }
  By defining own functions here, one can obtain an arbitrary degree of
  flexibility.
  
  Finally, \code{xaxis.labelFormat} is a \code{\link{strftime}}
  compatible formatting string., e.g. the default value is
  \code{"\%G\\n\\n\%OQ"}, which means ISO year and quarter (in roman
  letters) stacked on top of each other.
}

\value{
  \code{NULL} (invisibly). The function is called for its side effects.
}

\author{
  Michael H\enc{ö}{oe}hle with contributions by Sebastian Meyer
}

\seealso{
  the examples in \code{\link{stsplot_time1}} and \code{\link{plotHHH4_fitted1}}
}

\keyword{aplot}
