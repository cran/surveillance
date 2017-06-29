\name{meningo.age}
\alias{meningo.age}
\docType{data}
\title{Meningococcal infections in France 1985-1995}
\description{
  Monthly counts of meningococcal infections in France 1985-1995.
  Here, the data is split into 4 age groups (<1, 1-5, 5-20, >20).
}
\usage{data(meningo.age)}
\format{
  An object of class disProg with 156 observations in each one of 4 age groups.
  \describe{
    \item{week}{Number of month}
    \item{observed}{Matrix with number of counts in the corresponding month and age group}
    \item{state}{Boolean whether there was an outbreak -- dummy not implemented}
    \item{neighbourhood}{Neighbourhood matrix, all age groups are adjacent}
    \item{populationFrac}{Population fractions} 
  }
}
\source{
  ??
}
\examples{
data(meningo.age)
plot(meningo.age, title="Meningococcal infections in France 1985-95")
plot(meningo.age, as.one=FALSE)
}
\keyword{datasets}
