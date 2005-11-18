\name{plot.survRes}
\alias{plot.survRes}
\title{Plot Generation}
\description{Plotting of a survRes object.
}
\usage{
    plot.survRes(x, method="", disease="",
                        startyear = 2001, firstweek = 1, legend=TRUE, ...)
}
\arguments{
        \item{x}{Object of class \code{survRes}}
        \item{method}{Surveillance method to be used in title}
        \item{disease}{Name of disease in title}
        \item{startyear}{Year to begin the axis labeling (the year where
                the oldest data come from)}
        \item{firstweek}{Number of the first week of January in the first year (just for axis labeling reasons)}
        \item{legend}{Boolean indicating whether to add a legend}
        \item{...}{further arguments for the function \code{matplot}}
}
\value{
        \item{a plot}{showing the number of infected, the threshold for recognizing an
                outbreak, the computed alarm status and the defined alarm status.
        }
}
\author{M. Höhle, A. Riebler, C. Lang}
\examples{

    # Plotting of simulated data
    disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 400,
                                    A = 1, alpha = 1, beta = 0, phi = 0,
                                    frequency = 1, state = NULL, K = 5)
    # evaluate the timepoints defined by range using RKI 1
    control <- list(list(funcName = "rki1", range = 200:400))
    survResults <- algo.call(disProgObj, control = control)
    # plot the result
    plot(survResults[[1]], "RKI 1", "Simulation")
    plot(survResults[[1]], "RKI 1", "Simulation",
                                firstweek = 13, startyear = 2002)
    plot(survResults[[1]], "RKI 1", "Simulation", firstweek = 14)
}

\keyword{}


