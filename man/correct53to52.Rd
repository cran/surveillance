\name{correct53to52}
\alias{correct53to52}
\title{Data Correction from 53 to 52 weeks}
\description{Correction of data from 53 to 52 weeks a year
}
\usage{
    correct53to52(disProgObj, firstweek = 1)
}
\arguments{
        \item{disProgObj}{object of class disProg (including the observed and the state chain).}
        \item{firstweek}{the number of the first week in a year, default = 1 (if it starts
                                with the beginning of a year). Necessary, because the infected of week 53
                                and the infected of week 52 must be added.}
}
\value{
        \item{disProg}{a object \code{disProg} (disease progress) including a list of the observed
                 and the state chain (corrected to 52 weeks instead of 53 weeks a year)}
        }
}
\details{\code{\link{readData}} reads data with 53 weeks a year, but normally one year is said
        to have 52 weeks.
}
\seealso{\code{\link{readData}}
}
\examples{
    #This call correct53to52 automatically
    obj <- readData("k1",week53to52=TRUE)
    correct53to52(obj) # first entry is the first week of the year

    obj <- readData("n1",week53to52=FALSE)
    correct53to52(obj, firstweek = 5) # now it's assumed that the fifth
                                      # entry is the first week of the year
}

\keyword{}


