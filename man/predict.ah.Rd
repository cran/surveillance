\name{predict.ah}
\alias{predict.ah}
\title{Predictions from a HHH model}
\description{
  Use a \code{ah} object for prediction.
}
\usage{
  predict.ah(object,newdata, type=c("response","endemic"),...)
}
\arguments{
\item{object}{object of class \code{ah}} 
\item{newdata}{optionally, a disProgObject with which to predict; if omitted, the fitted mean is returned. }
\item{type}{the type of prediction required. The default is on the 
    scale of the response variable (endemic and epidemic part) 
    the alternative "endemic" returns only the endemic part 
    (i.e. \eqn{n_{it} * \nu_{it}})} 
\item{...}{not really used}
}
\value{
  matrix of values containing the mean \eqn{mu_{it}} for each region 
  and time point. 
}
\details{this function is still experimental}
\keyword{ models }



% algo.hhh.grid
