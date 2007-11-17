\name{toFileDisProg}
\alias{toFileDisProg}
\title{Writing of Disease Data}
\description{Writing of disease data (disProg object) into a file.
}
\usage{
    toFileDisProg(disProgObj, toFile)
}
\arguments{
        \item{disProgObj}{The disProgObj to save in file}
        \item{toFile}{The path and filename of the file to save}
}
\value{
        \item{file}{The file with the disease data}
}
\details{
  Writing of \code{disProg} object into a file as illustrated in the example.
}
\seealso{\code{\link{readData}}, \code{\link{sim.pointSource}}}
\examples{
    disProgObj <- sim.pointSource(length=200, K=1)
    toFileDisProg(disProgObj, "./simulation.txt")
    mydisProgObj <- readData("./simulation",sysPath=FALSE)
}

\keyword{file}


