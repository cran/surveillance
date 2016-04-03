\encoding{latin1}
\name{linelist2sts}
\alias{linelist2sts}

\title{
  Convert individual case information based on dates into an aggregated
  time series of counts
}

\description{
  The function is used to convert an individual line list of cases to an
  aggregated time series of counts based on event date information of the cases.
}

\usage{
linelist2sts(linelist,dateCol,
              aggregate.by=c("1 day", "1 week", "7 day", "1 week",
"1 month",  "3 month", "1 year"),
              dRange=NULL,
              epochInPeriodStr=switch(aggregate.by, "1 day"="1",
"1 week"="\%u", "1 month"="\%d","3 month"="\%q","1 year"="\%j"),
              startYearFormat=switch(aggregate.by,"1 day"="\%Y",
"7 day"="\%G", "1 week"="\%G","1 month"="\%Y","3 month"="\%Y","1 year"="\%Y"),
              startEpochFormat=switch(aggregate.by,"1 day"="\%j",
"7 day"="\%V", "1 week"="\%V", "1 month"="\%m", "3 month"="\%Q", "1 year"="1")
)
}
\arguments{
  \item{linelist}{
    A \code{data.frame} containing the line list of cases.
  }
  \item{dateCol}{A character string stating the column name in
    \code{linelist} which contains the event occurrence information (as a
    vector of \code{Date}s) which are to be temporally aggregated.
  }
  \item{aggregate.by}{Temporal aggregation level given as a string, see
  the \code{by} variable of the \code{\link{seq.Date}} function for
  further details.  
  }
  \item{dRange}{A vector containing the minimum and maximum date
    for doing the aggregation. If not specified these dates are
    extracted automatically by taking \code{range(D[,dateCol])} and
    adjust these according to \code{aggregate.by} (e.g. always first of
    a month).
  }
  \item{epochInPeriodStr}{\code{strptime} compatible format string to use for
    determining how a date is placed within the epoch. This is, e.g.,
    used to move the \code{dRange} epochs to the beginning of the
    period. Example: In case of weekly aggregation the "\%u" determines
    which day within the week (Monday is day 1) we have. See
    \code{\link{strptime}} for further details.
  }
  \item{startYearFormat}{\code{strptime} compatible format string to use for
    determining how the \code{start} entry of the \code{sts} object is
    generated. Usually the provided defaults are sufficient.}
  \item{startEpochFormat}{\code{strptime} compatible format string to use for
    determining how the \code{start} entry of the \code{sts} object is
    generated. Usually the provided defaults are sufficient.}
}

\details{
  The date range is automatically extended such that the starting and
  ending dates are always the first epoch within the period, i.e. for
  aggregation by week it is moved to Mondays. This is controlled by the
  \code{epochInPeriodStr} parameter.

  Please note that the formatting strings are implemented by the
  internal \code{formatDate} function, which uses the \code{strptime}
  formatting strings (ISO8601 are made Windows compatible) as well as
  formatting of quarters using "\%Q", "\%OQ" and "\%q".
}

\value{
  The function returns an object of class \code{"\linkS4class{sts}"}.
  The \code{freq} slot might not be appropriate.
}

\author{
Michael \enc{Höhle}{Hoehle}
}

\note{
  Some of the \code{strptime} ISO 8601 formatting strings have problems
  on windows. Hence, the surveillance internal \code{formatDate}
  function is used.
}

\seealso{
  See also \code{\link{seq.Date}}, \code{\link{strptime}}, \code{\link{isoWeekYear}}.
}

\examples{
#Load O104 outbreak data
data("husO104Hosp")

#Convert line list to an sts object
sts <- linelist2sts(husO104Hosp, dateCol="dHosp", aggregate.by="1 day")

#Check that the number of cases is correct
all.equal(sum(observed(sts)),nrow(husO104Hosp))

#Plot the result
plot(sts,xaxis.tickFreq=list("\%d"=atChange,"\%m"=atChange),
           xaxis.labelFreq=list("\%d"=at2ndChange),
           xaxis.labelFormat="\%d \%b",
           xlab="",las=2,cex.axis=0.8)
}

\keyword{models}
\keyword{optimize}
