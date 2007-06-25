\name{makePlot}
\alias{makePlot}
\title{Plot Generation}
\description{Just a test method.
}
\usage{
    makePlot(outputpath, data = "k1", method = "rki1",
                                name, disease, range = 157:339)
}
\arguments{
        \item{outputpath}{path for the storage}
        \item{data}{abbreviation of the disease-file}
        \item{method}{method to be called}
        \item{name}{name of the method}
        \item{disease}{disease name}
        \item{range}{range to plot}
}
\details{\code{makePlot} reads the data given in \code{data} using the function \code{readData},
        and the data are corrected to 52 weeks, enlarged using
        \code{enlargeData} and sendt to the surveillance system given in \code{method}.
        The system result is plotted and stored in \code{outputpath}.
}
\seealso{\code{\link{readData}}, \code{\link{correct53to52}}, \code{\link{enlargeData}},
                \code{\link{algo.call}}, \code{\link{plot.survRes}}
}
\author{M. Höhle, A. Riebler, C. Lang}
\examples{
    makePlot("./", "k1", "rki2", "RKI 2", "Kryptosporidiose")
}

\keyword{}

