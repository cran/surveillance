\name{create.disProg}
\alias{create.disProg}
\title{Creating an object of class disProg}
\description{
  Creates an object of class \code{disProg} from a vector with the weeknumber
  (week) and matrices with the observed number of counts (observed) and the
  respective state chains (state), where each column represents an individual
  time series. The matrices neighbourhood and populationFrac provide 
  information about neighbouring units and population proportions.
}
\usage{
  create.disProg(week, observed, state, start=c(2001,1), freq=52, 
  neighbourhood=NULL, populationFrac=NULL)
}
\arguments{
\item{week}{number of week}
\item{observed}{matrix with parallel time series of counts where rows are 
    time points and columns are the individual time series for unit/area
    \eqn{i, i=1,\ldots,m}}
\item{state}{matrix with corresponding states}
\item{start}{vector of length two denoting the year and the sample number (week, month, etc.) of the first observation}
\item{freq}{sampling frequency per year, i.e. 52 for weekly data, 12 for monthly data, 13 if 52 weeks are aggregated into 4 week blocks.}
\item{neighbourhood}{neighbourhood matrix \eqn{N} of dimension 
    \eqn{m \times m} with elements \eqn{n_{ij}=1} if units \eqn{i}
    and \eqn{j} are adjacent and 0 otherwise } 
\item{populationFrac}{matrix with corresponding population proportions} 
}
\value{\item{disProg}{object of class \code{disProg}  
  }
}
\details{}
\seealso{}
\author{M. Paul}
\examples{

  # create an univariate disProg object for the salmonella.agona data
  data(salmonella.agona)
  week <- nrow(salmonella.agona)
  salmonellaDisProg <- create.disProg(week=week, observed=salmonella.agona$observed,
  state=salmonella.agona$state)
  
  # plot salmonella cases
  title <- "Salmonella Agona cases in the UK"
  plot(salmonellaDisProg, title = title, xaxis.years=TRUE, legend=FALSE,
        startyear = 1990, firstweek = 1)
}
\keyword{}



% sumNeighbours
