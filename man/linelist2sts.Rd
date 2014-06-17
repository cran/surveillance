\encoding{latin1}
\name{linelist2sts}
\alias{linelist2sts}
%Internal functions

\title{
  Convert individual case information based on dates into an aggregated time series
}

\description{
  The function is used to convert an individual line list of cases to an
  aggregated time series based on date information of the cases.
}

\usage{
linelist2sts(linelist,dateCol,aggregate.by="1 week",dRange=NULL,
            startYearFormat=switch(aggregate.by,
"1 day"="\%V","7 day"="\%V","1 week"="\%V","1 month"="\%Y","3 month"="\%Y"),
           startEpochFormat=switch(aggregate.by,
"1 day"="\%j","7 day"="\%V","1 week"="\%V","1 month"="\%m","3 month"="\%Q")
)
}
\arguments{
  \item{linelist}{
    A \code{data.frame} containing the line list of cases.
  }
  \item{dateCol}{A character string stating the column name in
    \code{linelist} which contain the case data which are to be
    temporally aggregated.
  }
  \item{aggregate.by}{Temporal aggregation level given as a string, see
  the \code{by} variable of the \link{seq.Date} function for further details. 
  }
  \item{dRange}{A vector containing the minimum and maximum data to
    use. If not specified these dates are extracted automatically by
    taking \code{range(D[,dateCol])}.
  }
  \item{startYearFormat}{Strptime compatible format string to use for
    determining how the date string is generated. Usually the provided
  options are sufficient.}
  \item{startEpochFormat}{Strptime compatible format string to use for
    determining how the date string is generated. Usually the provided
  options are sufficient.}
}

\details{
  In case aggregation occurs by week the date range is automatically
  extended such that the starting and ending dates are mondays. This
  might not be an appropriate way in all situations this function is
  to be used.
}

\value{
  The function returns an object of class \code{"\linkS4class{sts}"}.
  The \code{freq} slot might not be appropriate.
}

\author{
Michael \enc{Höhle}{Hoehle}
}

\note{
  This implementation is still experimental, changes might occur in the future.
}

\seealso{
  See also \code{\link{seq.Date}}.
}

\examples{
#Load simulated outbreak data.
url <- paste("http://www.stat.uni-muenchen.de/~hoehle/",
             "teaching/moid2011/tutorials/cast-backnow/outbreak.txt",sep="")
D <- try(read.table(url,header=TRUE,colClasses=c("integer",rep("Date",3))))

if (!inherits(D, "try-error")) {
    #Convert line list to an sts object
    sts <- linelist2sts(D, dateCol="dOnset", aggregate.by="1 day")

    #Plot the result
    plot(sts,xaxis.tickFreq=list("\%d"=atChange,"\%m"=atChange),
             xaxis.labelFreq=list("\%d"=at2ndChange),
             xaxis.labelFormat="\%d \%b",legend.opts=NULL,
             xlab="",las=2,cex.axis=0.8)
}
}

\keyword{models}
\keyword{optimize}
