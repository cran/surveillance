\name{compMatrix.writeTable}
\alias{compMatrix.writeTable}
\encoding{latin1}

\title{Latex Table Generation}
\description{generates a latex table
}
\usage{
    compMatrix.writeTable(compMatrix)
}
\arguments{
        \item{compMatrix}{Matrix which includes quality values for every
                surveillance system.}
}
\value{
        \item{xtable}{
        Latex table of the entered matrix.
        }
}
\author{M. \enc{Höhle}{Hoehle}, A. Riebler, C. Lang}
\examples{
    ### First creates some tables ###

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
    control <- list( list(funcName = "rki1", range = range),
                    list(funcName = "rki2", range = range),
                    list(funcName = "rki3", range = range)
                )

    ### This are single compMatrices
    compMatrix1 <- algo.compare(algo.call(disProgObj1, control=control))
    compMatrix2 <- algo.compare(algo.call(disProgObj2, control=control))
    compMatrix3 <- algo.compare(algo.call(disProgObj3, control=control))

    ### This is a summary compMatrix
    sumCompMatrix <- algo.summary( list(a=compMatrix1,
                                            b=compMatrix2, c=compMatrix3) )

    ### Now show the latextable from the single compMatrix compMatrix1
    compMatrix.writeTable(compMatrix1)

    ### Now show the latextable from the summary compMatrix
    compMatrix.writeTable(sumCompMatrix)
}

\keyword{print}

