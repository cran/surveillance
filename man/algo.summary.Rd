\name{algo.summary}
\alias{algo.summary}
\title{Summary Table Generation for Several Disease Chains}
\description{
  Summary table generation for several disease chains.
}
\usage{
algo.summary(compMatrices)
}
\arguments{
  \item{compMatrices}{list of matrices constructed by algo.compare.}
}
\value{
  a matrix summing up the singular input matrices
}
\details{
  As lag the mean of all single lags is returned. TP values, FN values,
  TN values and FP values are summed up. \code{dist}, \code{sens} and
  \code{spec} are new computed on the basis of the new TP value, FN value,
  TN value and FP value.
}
\seealso{\code{\link{algo.compare}}, \code{\link{algo.quality}}}
\examples{
# Create a test object
disProgObj1 <- sim.pointSource(p = 0.99, r = 0.5, length = 400,
                               A = 1, alpha = 1, beta = 0, phi = 0,
                               frequency = 1, state = NULL, K = 1.7)
disProgObj2 <- sim.pointSource(p = 0.99, r = 0.5, length = 400,
                               A = 1, alpha = 1, beta = 0, phi = 0,
                               frequency = 1, state = NULL, K = 5)
disProgObj3 <- sim.pointSource(p = 0.99, r = 0.5, length = 400,
                               A = 1, alpha = 1, beta = 0, phi = 0,
                               frequency = 1, state = NULL, K = 17)

# Let this object be tested from any methods in range = 200:400
range <- 200:400
control <- list(list(funcName = "rki1", range = range),
                list(funcName = "rki2", range = range),
                list(funcName = "rki3", range = range))

compMatrix1 <- algo.compare(algo.call(disProgObj1, control=control))
compMatrix2 <- algo.compare(algo.call(disProgObj2, control=control))
compMatrix3 <- algo.compare(algo.call(disProgObj3, control=control))

algo.summary( list(a=compMatrix1, b=compMatrix2, c=compMatrix3) )
}
\keyword{print}
