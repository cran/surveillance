\name{enlargeData}
\alias{enlargeData}
\title{Data Enlargement}
\description{Enlargement of data which is too short for a surveillance method to evaluate.
}
\usage{
    enlargeData(disProgObj, range = 1:156, times = 1)
}
\arguments{
\item{disProgObj}{object of class disProg (including the observed and the state chain).}
\item{range}{range of already existing data (\code{state}, \code{observed})
    which should be used for enlargement.}
\item{times}{number of times to enlarge.}
}
\value{
\item{disProg}{a object \code{disProg} (disease progress) including a list of the observed
    and the state chain (extended with cyclic data generation)}
}
\details{
        \code{observed} and \code{state} are enlarged in the way that the part \code{range}
        of \code{observed} and \code{state}
        is repeated \code{times} times in front of \code{observed} and \code{state}.
        Sometimes it's useful to care for the cyclic property of the timeseries, so as default
        we enlarge observed and state once with the first three existing years, assuming a year
        has 52 weeks.
}
\seealso{\code{\link{readData}}
}
\examples{
    obj <- readData("k1")

    enlargeData(obj) # enlarge once with part 1:156
    enlargeData(obj, 33:36, 10) # enlarge 10 times with part 33:36
}

\keyword{utilities}

