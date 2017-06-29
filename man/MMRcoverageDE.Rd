\name{MMRcoverageDE}
\alias{MMRcoverageDE}
\docType{data}
\title{MMR coverage levels in the 16 states of Germany}
\description{
  Coverage levels at school entry for the first and second dose
  of the combined measles-mumps-rubella (MMR) vaccine in 2006, estimated 
  from children presenting vaccination documents at school entry examinations.
}
\usage{data(MMRcoverageDE)}
\format{
  A \code{data.frame} containing 19 rows and 5 columns with variables
  \describe{
    \item{state}{Names of states: the 16 federal states are
    followed by the total of Germany, as well as the total 
    of West and East Germany.}
    \item{nOfexaminedChildren}{Number of children examined.}
    \item{withVaccDocument}{Percentage of children who presented vaccination documents.}
    \item{MMR1}{Percentage of children with vaccination documents, who
      received at least 1 dose of MMR vaccine.}
    \item{MMR2}{Percentage of children with vaccination documents, who
      received at least 2 doses of MMR vaccine.}
  } 
  
  Coverage levels were derived from vaccination documents presented
  at medical examinations, which are conducted by local health authorities
  at school entry each year. Records include information about the receipt 
  of 1st and 2nd doses of MMR, but no information about dates.
  Note that information from children who did not present a vaccination 
  document on the day of the medical examination, is not included in
  the estimated coverage. 
}
\source{
  Robert Koch-Institut (2008) Zu den Impfquoten bei den 
    Schuleingangsuntersuchungen in Deutschland 2006. 
    Epidemiologisches Bulletin, \bold{7}, 55-57
}
\seealso{\code{\link{measlesDE}}}
\references{
  Herzog, S.A., Paul, M. and Held, L. (2011) Heterogeneity in vaccination
    coverage explains the size and occurrence of measles epidemics in
    German surveillance data. Epidemiology and Infection, \bold{139}, 505--515.
}
\keyword{datasets}
