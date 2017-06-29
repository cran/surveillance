\name{residuals.ah}
\alias{residuals.ah}
\alias{residuals.ahg}
\title{Residuals from a HHH model}
\description{
  Extracts model residuals from a \code{ah} or \code{ahg} object.
}
\usage{
\method{residuals}{ah}(object, type=c("deviance","pearson"), \dots)
}
\arguments{
\item{object}{object of class \code{ah} or \code{ahg} } 
\item{type}{the type of residuals which should be returned. 
  The alternatives are "deviance" (default) and "pearson"} 
\item{\dots}{not really used}
}
\value{
  matrix with residuals for each region 
  and time point. 
}
\note{This function is experimental!}
\keyword{models}
