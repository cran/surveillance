\name{salmonella.agona}
\alias{salmonella.agona}
\docType{data}
\title{Salmonella Agona cases in the UK 1990-1995}
\description{
  Reported number of cases of the Salmonella Agona serovar in the UK
  1990-1995. Note however
  that the counts do not correspond exactly to the ones used by
  Farrington et. al (1996).
}
\usage{data(salmonella.agona)}
\format{
  A data frame with 312 observations on the following 2 variables.
  \describe{
    \item{week}{First four digits are the year, last two the week
      number within that year}
    \item{observed}{Number of counts in the corresponding week}
    \item{state}{Boolean whether there was an outbreak -- dummy not implemented.}
  }
}
\source{A statistical algorithm for the early detection of outbreaks of infectious disease, Farrington, C.P., Andrews, N.J, Beale A.D. and Catchpole, M.A. (1996). , J. R. Statist. Soc. A, 159, 547-563.}
\examples{
data(salmonella.agona)
plot(salmonella.agona$observed,type="l",ylab="counts",xlab="")

}
\keyword{datasets}

