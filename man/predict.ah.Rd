\name{predict.ah}
\alias{predict.ah}
\alias{predict.ahg}
\encoding{latin1}
\title{Predictions from a HHH model}
\description{
  Use a \code{ah} or \code{ahg} object for prediction.
}
\usage{
\method{predict}{ah}(object,newdata=NULL, 
   type=c("response","endemic","epi.own","epi.neighbours"), \dots)
}
\arguments{
\item{object}{object of class \code{ah} or \code{ahg} } 
\item{newdata}{optionally, a disProgObject with which to predict; if omitted, the fitted mean is returned. }
\item{type}{the type of prediction required. The default is on the 
    scale of the response variable (endemic and epidemic part). 
    The alternative "endemic" returns only the endemic part 
    (i.e. \eqn{n_{it} \nu_{it}}{n_it * \nu_it}), "epi.own" and "epi.neighbours"
    return the epidemic part (i.e. \eqn{\lambda_i y_{i,t}}{\lambda_i * y_i,t}
    and \eqn{\phi_i \sum_{j \sim i} y_{j,t-1}}{\phi_i * \sum_(j ~ i) y_j,t-1} )} 
\item{...}{not really used}
}
\value{
  matrix of values containing the mean \eqn{\mu_{it}}{\mu_it} for each region 
  and time point. 
}
\details{this function is still experimental}
\keyword{ models }



% residuals.ah
