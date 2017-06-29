\name{algo.quality}
\alias{algo.quality}
\title{Computation of Quality Values for a Surveillance System Result}
\description{
  Computation of the quality values for a surveillance system output.
}
\usage{
algo.quality(sts, penalty = 20)
}
\arguments{
\item{sts}{object of class \code{survRes} or \code{sts}, which includes the state chain and
    the computed alarm chain}
\item{penalty}{the maximal penalty for the lag}
}
\value{
  a list of quality values:
  \item{TP}{Number of correct found outbreaks.}
  \item{FP}{Number of false found outbreaks.}
  \item{TN}{Number of correct found non outbreaks.}
  \item{FN}{Number of false found non outbreaks.}
  \item{sens}{True positive rate, meaning TP/(FN + TP).}
  \item{spec}{True negative rate, meaning TN/(TN + FP).}
  \item{dist}{Euclidean distance between (1-spec, sens) to (0,1).}
  \item{lag}{Lag of the outbreak recognizing by the system.}
}
\details{
  The lag is defined as follows:
  In the state chain just the beginnings of an outbreak chain (outbreaks directly
  following each other) are considered. In the alarm chain, the range from the beginning
  of an outbreak until \eqn{min(next outbreak beginning,\code{penalty})} timepoints is considered. The \code{penalty} timepoints were
  chosen, to provide an upper bound on the penalty for not discovering an outbreak. Now the difference between the first alarm by the system and the defined beginning is denoted ``the lag'' Additionally outbreaks found by the system are not
  punished. At the end, the mean of the lags for every outbreak chain is returned
  as summary lag.
}
\seealso{\code{\link{algo.compare}}}
\examples{
# Create a test object
disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 200, A = 1,
                              alpha = 1, beta = 0, phi = 0,
                              frequency = 1, state = NULL, K = 1.7)

# Let this object be tested from rki1
survResObj <- algo.rki1(disProgObj, control = list(range = 50:200))

# Compute the quality values
algo.quality(survResObj)
}
\keyword{misc}
