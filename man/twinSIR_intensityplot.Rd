\encoding{latin1}
\name{twinSIR_intensityplot}
\alias{plot.twinSIR}
\alias{intensityplot.twinSIR}
\alias{intensityplot.simEpidata}

\title{
  Plotting Paths of Infection Intensities for \code{twinSIR} Models
}

\description{
\code{\link{intensityplot}} methods to plot the evolution of the total infection
intensity, its epidemic proportion or its endemic proportion over time.
The default \code{plot} method for objects of class \code{"twinSIR"}
is just a wrapper for the \code{intensityplot} method.
The implementation is illustrated in Meyer et al. (2017, Section 4),
see \code{vignette("twinSIR")}.
}

\usage{
\method{plot}{twinSIR}(x, which = c("epidemic proportion", "endemic proportion",
     "total intensity"), ...)

\method{intensityplot}{twinSIR}(x, which = c("epidemic proportion", "endemic proportion",
              "total intensity"), aggregate = TRUE, theta = NULL,
              plot = TRUE, add = FALSE, rug.opts = list(), ...)

\method{intensityplot}{simEpidata}(x, which = c("epidemic proportion", "endemic proportion",
              "total intensity"), aggregate = TRUE, theta = NULL,
              plot = TRUE, add = FALSE, rug.opts = list(), ...)
}

\arguments{
  \item{x}{
    an object of class \code{"\link{twinSIR}"} (fitted model) or
    \code{"\link{simEpidata}"} (simulated \code{twinSIR} epidemic),
    respectively. 
  }
  \item{which}{
    \code{"epidemic proportion"}, \code{"endemic proportion"},
    or \code{"total intensity"}.  Partial matching is applied.  Determines
    whether to plot the path of the total intensity \eqn{\lambda(t)} or its 
    epidemic or endemic proportions
    \eqn{\frac{e(t)}{\lambda(t)}}{e(t)/lambda(t)} or
    \eqn{\frac{h(t)}{\lambda(t)}}{h(t)/lambda(t)}. 
  }
  \item{aggregate}{
    logical. Determines whether lines for all individual infection
    intensities should be drawn (\code{FALSE}) or their sum only
    (\code{TRUE}, the default).
  }
  \item{theta}{
    numeric vector of model coefficients.  If \code{x} is of class
    \code{"twinSIR"}, then \code{theta = c(alpha, beta)}, where \code{beta}
    consists of the coefficients of the piecewise constant log-baseline function
    and the coefficients of the endemic (\code{cox}) predictor.  If \code{x} is
    of class \code{"simEpidata"}, then \code{theta = c(alpha, 1, betarest)}, 
    where 1 refers to the (true) log-baseline used in the simulation and 
    \code{betarest} is the vector of the remaining coefficients of the endemic 
    (\code{cox}) predictor.
    The default (\code{NULL}) means that the fitted or true parameters,
    respectively, will be used.
  }
  \item{plot}{
    logical indicating if a plot is desired, defaults to \code{TRUE}.
    Otherwise, only the data of the plot will be returned.  Especially with
    \code{aggregate = FALSE} and many individuals one might e.g. consider to 
    plot a subset of the individual intensity paths only or do some further
    calculations/analysis of the infection intensities.
  }
  \item{add}{
    logical.  If \code{TRUE}, paths are added to the current plot, using
    \code{lines}.
  }
  \item{rug.opts}{
    either a list of arguments passed to the function \code{\link{rug}}, or
    \code{NULL} (or \code{NA}), in which case no \code{rug} will be plotted.
    By default, the argument \code{ticksize} is set to 0.02 and \code{quiet}
    is set to \code{TRUE}.  Note that the argument \code{x} of the
    \code{rug()} function, which contains the
    locations for the \code{rug} is fixed internally and can not be modified.
    The locations of the rug are the time points of infections.
  }
  \item{\dots}{
    For the \code{plot.twinSIR} method, arguments passed to
    \code{intensityplot.twinSIR}.
    For the \code{intensityplot} methods, further graphical parameters
    passed to the function \code{\link{matplot}}, 
    e.g. \code{lty}, \code{lwd}, \code{col}, \code{xlab}, \code{ylab} and
    \code{main}.  Note that the \code{matplot} arguments \code{x}, \code{y},
    \code{type} and \code{add} are implicit and can not be specified here.
  }
}

\value{
  numeric matrix with the first column \code{"stop"} and as many rows as there
  are \code{"stop"} time points in the event history \code{x}.  The other 
  columns depend on the argument \code{aggregate}: if \code{TRUE}, there 
  is only one other column named \code{which}, which contains the values of 
  \code{which} at the respective \code{"stop"} time points.  Otherwise, if
  \code{aggregate = FALSE}, there is one column for each individual, each of
  them containing the individual \code{which} at the respective \code{"stop"}
  time points.
}

\references{
  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
}

\author{
Sebastian Meyer
}

\seealso{
\code{\link{twinSIR}} for a description of the intensity model, and
\code{\link{simulate.twinSIR}} for the simulation of epidemic data
according to a \code{twinSIR} specification.
}

\examples{
data("hagelloch")
plot(hagelloch)

# a simplistic twinSIR model
fit <- twinSIR(~ household, data = hagelloch)

# overall total intensity
plot(fit, which = "total")

# overall epidemic proportion
epi <- plot(fit, which = "epidemic", ylim = c(0, 1))
head(epi)
# add overall endemic proportion = 1 - epidemic proportion
ende <- plot(fit, which = "endemic", add = TRUE, col = 2)
legend("topleft", legend = "endemic proportion", lty = 1, col = 2, bty = "n")

# individual intensities
tmp <- plot(fit, which = "total", aggregate = FALSE,
    col = rgb(0, 0, 0, alpha = 0.1),
    main = expression("Individual infection intensities " *
        lambda[i](t) == Y[i](t) \%.\% (e[i](t) + h[i](t))))
# return value: matrix of individual intensity paths
str(tmp)

# plot intensity path only for individuals 3 and 99
matplot(x = tmp[,1], y = tmp[,1+c(3,99)], type = "S",
        ylab = "Force of infection", xlab = "time",
        main = expression("Paths of the infection intensities " *
                          lambda[3](t) * " and " * lambda[99](t)))
legend("topright", legend = paste("Individual", c(3,99)),
       col = 1:2, lty = 1:2)
}

\keyword{hplot}
\keyword{aplot}
\keyword{dplot}
\keyword{methods}
