\name{magic.dim}
\alias{magic.dim}
\title{Compute Suitable k1 x k2 Layout for Plotting}
\description{
       For a given number \code{k}, \code{magic.dim} provides a vector 
       containing two elements, the number of rows (k1) and columns (k2),
       respectively, which can be used to set the 
       dimension of a single graphic device so that k1*k2 plots can be 
       drawn by row (or by column) on the device.
}
\usage{
      magic.dim(k)
}
\arguments{
      \item{k}{an integer}
}
\value{numeric vector with two elements}
\seealso{
\code{\link{primeFactors}} and \code{\link{bestCombination}} which are
internally used to complete the task.

\code{\link{n2mfrow}} is a similar function from package \pkg{grDevices}.
}
\keyword{dplot}
\keyword{utilities}
