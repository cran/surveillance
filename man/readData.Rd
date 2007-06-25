\name{readData}
\alias{readData}
\title{Reading of Disease Data}
\description{Reading of disease data. In the package disease data are saved in a file \code{<abb>.txt} containing three columns -- the weeknumber (week), the observed number of counts (observed) and a state (state). The data are read using \code{read.table(...,header=T)}, hence the file has to contain a header.

}
\usage{
    readData(abb,week53to52=TRUE,sysPath=TRUE)
}
\arguments{
        \item{abb}{abbreviation of the diseasename.}
        \item{week53to52}{Boolean indicating whether to convert RKI 53 Weeks System to 52 weeks a year} 
          \item{sysPath}{Boolean, if \code{TRUE} then R automatically looks in the data directory of the surveillance package.}
}
\value{
        \item{disProg}{a object \code{disProg} (disease progress) including a list of the observed and the state chain.}
}
\details{
  This function is only kept for backwards compability. As of 0.9-2 all data should be read with \code{data}.
}
\seealso{\code{\link{m1}}, \code{\link{m2}}, \code{\link{m3}}, \code{\link{m4}}, \code{\link{m5}},
        \code{\link{q1_nrwh}}, \code{\link{q2}}, \code{\link{s1}}, \code{\link{s2}}, \code{\link{s3}},
        \code{\link{k1}}, \code{\link{n1}}, \code{\link{n2}}, \code{\link{h1_nrwrp}}
}
\examples{
    readData("m5")

    #To bring a single vector of counts into a format, which can be
    #handled by readData. Assume ``counts'' is a vector of counts.
    counts <- rpois(100,20)
    counts <- data.frame("week"=1:length(counts),"observed"=counts,"state"=rep(0,length(counts)))
    write(c("week","observed","state"),file="disease.txt",ncol=3)
    write(t(as.matrix(counts)),file="disease.txt",ncol=3,append=TRUE)
    disease <- readData("disease",week53to52=FALSE,sysPath=FALSE)

}

\keyword{}


