\name{deleval}
\alias{deleval}
\docType{data}
\title{Surgical Failures Data}
\description{
  The dataset from Steiner et al. (1999) on A synthetic dataset from the Danish meat inspection -- useful
  for illustrating the beta-binomial CUSUM.
}
\usage{data(abattoir)}
\details{
Steiner et al. (1999) use data from de Leval et al. (1994) to
illustrate monitoring of failure rates of a surgical procedure
for a bivariate outcome. 

Over a period of six years an arterial switch operation was performed
on 104 newborn babies. Since the death rate from this surgery was
relatively low the idea of surgical "near miss" was introduced. It is
defined as the need to reinstitute cardiopulmonary bypass after a trial
period of weaning. The object of class \code{sts} contains the
recordings of near misses and deaths from the surgery for the 104
newborn babies of the study.

The data could also be handled by a multinomial CUSUM model.
}
\seealso{\code{\link{pairedbinCUSUM}}}
\examples{
data("deleval")
plot(deleval, xaxis.labelFormat=NULL,ylab="Response",xlab="Patient number")
}
\references{
  Steiner, S. H., Cook, R. J., and Farewell, V. T. (1999), Monitoring
  paired binary surgical outcomes using cumulative sum charts,
  Statistics in Medicine, 18, pp. 69--86.
  
  De Leval, Marc R., Franiois, K., Bull, C., Brawn, W. B. and
  Spiegelhalter, D. (1994),  Analysis of a cluster of
  surgical failures, Journal of Thoracic and Cardiovascular Surgery,
  March, pp. 914--924.
}
\keyword{datasets}
