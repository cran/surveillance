\name{epidata_plot}
\alias{plot.epidata}
\alias{plot.summary.epidata}
\alias{stateplot}

\title{
Plotting the Evolution of an Epidemic
}

\description{
  Functions for plotting the evolution of epidemics. The \code{\link{plot}}
  methods for \code{\link{class}}es \code{"\link{epidata}"} and
  \code{"summary.epidata"} plots the numbers of susceptible, infectious and
  recovered (= removed) individuals by step functions along the time axis.  The
  function \code{stateplot} shows individual state changes along the time axis.
}

\usage{
\method{plot}{summary.epidata}(x,
     lty = c(2, 1, 3), lwd = 2,
     col = c("#1B9E77", "#D95F02", "#7570B3"), col.hor = col, col.vert = col,
     xlab = "Time", ylab = "Number of individuals",
     xlim = NULL, ylim = NULL, legend.opts = list(), do.axis4 = NULL,
     panel.first = grid(), rug.opts = list(),
     which.rug = c("infections", "removals", "susceptibility", "all"), ...)
\method{plot}{epidata}(x, ...)

stateplot(x, id, ...)
}

\arguments{
  \item{x}{
    an object inheriting from class \code{"epidata"} or
    \code{"summary.epidata"}.  In the former case, its summary is calculated
    and the function continues as in the latter case.  The \code{plot} method
    for class \code{"epidata"} is a simple wrapper for
    \code{plot.summary.epidata} implemented as \code{plot(summary(x, ...))}.
  }
  \item{lty, lwd}{
    vectors of length 3 containing the line types and widths, respectively, for
    the numbers of susceptible, infectious and removed individuals (in this 
    order).  By default, all lines have width 1 and the line types are dashed
    (susceptible), solid (infectious) and dotted (removed), respectively.  To 
    omit the drawing of a specific line, just set the corresponding entry in 
    \code{lty} to 0.  The vectors are recycled if necessary.  For information 
    about the different \code{lty} and \code{lwd} codes, see the help pages 
    of \code{\link{par}}.
  }
  \item{col, col.hor, col.vert}{
    vectors of length 3 containing the line colors for the numbers of
    susceptible, infectious and removed individuals (in this order). 
    \code{col.hor} defines the color for the horizontal parts of the step
    function, whilst \code{col.vert} defines the color for its vertical parts.
    The argument \code{col} is just short for \code{col.hor = col} and
    \code{col.vert = col}.  The default \code{col} vector corresponds to
    \code{\link[RColorBrewer]{brewer.pal}("Dark2",n=3)} from the
    \pkg{RColorBrewer} package.
    The vectors are recycled if necessary.  For information about the possible
    values of \code{col}, see the help pages of \code{\link{par}}.
  }
  \item{xlab, ylab}{
    axis labels, default to "Time" and "Number of individuals", respectively.
  }
  \item{xlim, ylim}{
    the x and y limits of the plot in the form \code{c(xmin, xmax)} and
    \code{c(ymin, ymax)}, respectively.  By default, these are chosen adequately
    to fit the time range of the epidemic and the number of individuals.
  }
  \item{legend.opts}{
    if this is a list (of arguments for the \code{\link{legend}} function),
    a legend will be plotted.  The defaults are as follows:
    \describe{
      \item{\code{x}:}{\code{"topright"}}
      \item{\code{inset}:}{\code{c(0,0.02)}}
      \item{\code{legend}:}{\code{c("susceptible", "infectious", "removed")}}
      \item{\code{lty},\code{lwd},\code{col}:}{same as the arguments
        \code{lty}, \code{lwd}, and \code{col.hor} of the main function}
      \item{\code{bty}:}{\code{"n"}}
    }
  }
  \item{do.axis4}{
    logical indicating if the final numbers of susceptible and removed
    individuals should be indicated on the right axis.  The default \code{NULL}
    means \code{TRUE}, if \code{x} represents a SIR epidemic and \code{FALSE}
    otherwise, i.e. if the epidemic is SI, SIS or SIRS.
  }
  \item{panel.first}{
    an expression to be evaluated after the plot axes are set up but before
    any plotting takes place.  By default, a standard grid is drawn.
  }
  \item{rug.opts}{
    either a list of arguments passed to the function \code{\link{rug}} or
    \code{NULL} (or \code{NA}), in which case no \code{rug} will be plotted.
    By default, the argument \code{ticksize} is set to 0.02, \code{col}
    is set to the color according to \code{which.rug} (black if this is
    \code{"all"}), and \code{quiet} is set to \code{TRUE}.
    Note that the argument \code{x}, which contains the
    locations for the \code{rug} is fixed internally and can not be modified.
    The argument \code{which.rug} (see below) determines the locations to mark.
  }
  \item{which.rug}{
    By default, tick marks are drawn at the time points of infections.  
    Alternatively, one can choose to mark only \code{"removals"},
    \code{"susceptibilities"} (i.e. state change from R to S) or
    \code{"all"} events.
  }
  \item{id}{
    single character string or factor of length 1 specifying the individual for
    which the \code{stateplot} should be established.
  }
  \item{\dots}{
    For \code{plot.summary.epidata}: further graphical parameters passed to 
    \code{plot}, \code{lines} and \code{axis}, e.g. \code{main}, \code{las}, 
    \code{cex.axis} (etc.) and \code{mgp}.\cr
    For \code{plot.epidata}: arguments passed to \code{plot.summary.epidata}.\cr
    For \code{stateplot}: arguments passed to \code{\link{plot.stepfun}} or
    \code{\link{plot.function}} (if \code{id} had no events during the
    observation period).  By default, \code{xlab="time"}, \code{ylab="state"},
    \code{xlim=attr(x,"timeRange")}, \code{xaxs="i"} and \code{do.points=FALSE}.
  }
}

\value{
  \code{plot.summary.epidata} (and \code{plot.epidata}) invisibly returns the
  matrix used for plotting, which contains the evolution of the three
  counters.\cr
  \code{stateplot} invisibly returns the function, which was plotted,
  typically of class \code{"stepfun"}, but maybe of class \code{"function"},
  if no events have been observed for the individual in question (then the
  function always returns the initial state).  The vertical axis of
  \code{stateplot} can range from 1 to 3, where 1 corresponds to
  \emph{S}usceptible, 2 to \emph{I}nfectious and 3 to \emph{R}emoved.
}

\author{
Sebastian Meyer
}

\seealso{
\code{\link{summary.epidata}} for the data, on which the plots are based.
\code{\link{animate.epidata}} for the animation of epidemics.
}

\examples{
data("hagelloch")
(s <- summary(hagelloch))

# rudimentary stateplot
stateplot(s, id = "187")

# evolution of the epidemic
plot(s)
}

\keyword{hplot}
\keyword{methods}
\keyword{spatial}
