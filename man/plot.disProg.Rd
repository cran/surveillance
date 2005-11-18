\name{plot.disProg}
\alias{plot.disProg}
\title{Plot Generation of the Observed and the defined Outbreak State of a Timeseries}
\description{Plotting of a disProg object.
}
\usage{
    plot.disProg(x, title = "", startyear = 2001, firstweek = 1, ...)
}
\arguments{
        \item{x}{Object of class \code{disProg}}
        \item{title}{Plot title}
        \item{startyear}{Year to begin the axis labeling (the year where
                the oldest data come from)}
        \item{firstweek}{Number of the first week of January in the first year (just for axis labeling grounds)}
        \item{...}{further arguments for the function \code{matplot}}
}
\value{
        \item{a plot}{showing the number of infected and the defined alarm status
                for a timeseries created by simulation or given in data.
        }
}
\author{M. Höhle, A. Riebler, C. Lang}
\examples{

    # Plotting of simulated data
    disProgObj <- sim.pointSource(p = 0.99, r = 0.5, length = 208,
                                    A = 1, alpha = 1, beta = 0, phi = 0,
                                    frequency = 1, state = NULL, K = 5)
    # plot the simulated disease with the defined outbreaks
    plot(disProgObj)
    title <- "Number of Infected and Defined Outbreak Positions for Simulated Data"
    plot(disProgObj, title = title)
    plot(disProgObj, title = title,
            startyear = 1999, firstweek = 13)
    plot(disProgObj, title = title,
            startyear = 1999, firstweek = 14)
}

\keyword{}

