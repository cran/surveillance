\name{toLatex.sts}
\docType{methods}
\alias{toLatex.sts}
\alias{toLatex,sts-method}

\title{\code{toLatex}-Method for \code{"sts"} Objects}

\description{
Convert \code{"\linkS4class{sts}"} objects to a
character vector with LaTeX markup.
}

\usage{
\S4method{toLatex}{sts}(object, caption = "",label=" ", columnLabels = NULL,
        subset = NULL, 
        alarmPrefix = "\\\\textbf{\\\\textcolor{red}{",
        alarmSuffix = "}}", ubColumnLabel = "UB", ...)
}

\arguments{
\item{object}{an \code{"\linkS4class{sts}"} object.}
\item{caption}{A caption for the table. Default is the empty string.}
\item{label}{A label for the table. Default is the empty string.}
\item{columnLabels}{A list of labels for each column of the resulting table. Default is NULL}
\item{subset}{A range of values which should be displayed. If Null, then all data in the sts objects will be displayed. Else only a subset of 
data. Therefore range needs to be a numerical vector of indexes from 1 to length(@observed).}
\item{alarmPrefix}{A latex compatible prefix string wrapped around a table cell iff there is an alarm;i.e. alarm = TRUE}
\item{alarmSuffix}{A latex compatible suffix string wrapped around a table cell iff there is an alarm;i.e. alarm[i,j] = TRUE}
\item{ubColumnLabel}{The label of the upper bound column; default is
  \"UB\".}
\item{\dots}{further arguments passed to \code{\link{print.xtable}}.}
}

\value{
 An object of class \code{\link[=toLatex]{"Latex"}}.
}

\examples{
# Create a test object
data("salmonella.agona")

# Create the corresponding sts object from the old disProg object
salm <- disProg2sts(salmonella.agona)

control <- list(range=(260:312),
                noPeriods=1,populationOffset=FALSE,
                fitFun="algo.farrington.fitGLM.flexible",
                b=4,w=3,weightsThreshold=1,
                pastWeeksNotIncluded=3,
                pThresholdTrend=0.05,trend=TRUE,
                thresholdMethod="delta",alpha=0.1)
salm <- farringtonFlexible(salm,control=control)
toLatex(salm)
}

\author{Dirk Schumacher}
\keyword{print}
