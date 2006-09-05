\name{plot.disProg}
\alias{plot.disProg}
\alias{plot.disProg.one}
\title{Plot Generation of the Observed and the defined Outbreak States of a 
         (multivariate) time series}
\description{Plotting of a disProg object.
}
\usage{
    plot.disProg(x, title = "", xaxis.years=TRUE, startyear = 2001, 
         firstweek = 1, as.one=TRUE, same.scale=TRUE, legend=TRUE, ...)
    plot.disProg.one(x, title = "", xaxis.years=TRUE, startyear = 2001, firstweek = 1, ylim=NULL, legend=TRUE, xlegpos = 1/4, ylegpos = 1,...)
}
\arguments{
        \item{x}{Object of class \code{disProg}}
        \item{title}{Plot title}
        \item{xaxis.years}{if \code{true}, the x axis is labeled}
        \item{startyear}{Year to begin the axis labeling (the year where
                the oldest data come from)}
        \item{firstweek}{Number of the first week of January in the first year 
                (just for axis labeling grounds)}
        \item{as.one}{if \code{true} all individual time series are shown in 
              one plot}
        \item{same.scale}{if \code{true} all plots have same scale}
        \item{ylim}{range of y axis}
        \item{legend}{if \code{true} a legend is shown}
        \item{xlegpos}{x position of the legend as fraction of \code{max(xlim)}}
        \item{ylegpos}{y position of legend as fraction of \code{max(ylim)}}
        \item{...}{further arguments for the function \code{matplot}}
}
\value{
        \item{a plot}{showing the number of infected and the defined alarm status
                for a time series created by simulation or given in data 
                either in one single plot or in several plots for each individual
                time series. }
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
    plot(disProgObj, title = title, xaxis.years=TRUE,
            startyear = 1999, firstweek = 13)
    plot(disProgObj, title = title, xaxis.years=TRUE,
            startyear = 1999, firstweek = 14)
    
    # Plotting of measels data
    data(measels.weser)
    # one plot
    plot(measels.weser, title = "measels cases in the district Weser-Ems",
            xaxis.years=TRUE, startyear= 2001, firstweek=1)
    # plot cases for each "Kreis" 
    plot(measels.weser, same.scale=TRUE, as.one=FALSE)
}

\keyword{}




