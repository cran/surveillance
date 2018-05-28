\encoding{latin1}
\name{stsplot_spacetime}
\alias{stsplot_spacetime}
\title{
  Map of Disease Incidence
}

\description{
  For each period (row) or for the overall period of the
  \code{observed} matrix of the \code{"\linkS4class{sts}"} object, a map
  showing the counts by region is produced.
  It is possible to redirect the output into files, e.g., to generate an
  animated GIF.
}

\usage{
stsplot_spacetime(x, type, legend = NULL, opts.col = NULL, labels = TRUE,
                  wait.ms = 250, cex.lab = 0.7, verbose = FALSE,
                  dev.printer = NULL, ...)
}

\arguments{
  \item{x}{
    an object of class \code{"\linkS4class{sts}"}.
  }
  \item{type}{
    a formula (see \code{\link{stsplot}}). For a map aggregated over
    time (no animation), use \code{observed ~ 1 | unit}, otherwise
    \code{observed ~ 1 | unit * time}.
  }
  \item{legend}{
    An object of type \code{list} containing the following items used
    for coloring
    \itemize{
      \item{dx}{position increments in x direction}
      \item{dy}{position increments in y direction}
      \item{x}{position in x}
      \item{y}{position in y}
      \item{once}{\code{Boolean} - if \code{TRUE} then only shown once}
    }
    If \code{NULL} then a default legend is used.
  }
  \item{opts.col}{
    A list containing the two elements
    \itemize{
      \item{ncolors}{Number of colors to use for plotting}
      \item{use.color}{\code{Boolean} if \code{TRUE} then colors
	will be used in the palette, otherwise grayscale}
    }
  }
  \item{labels}{\code{Boolean} whether to add labels
  }
  \item{wait.ms}{Number of milliseconds to wait between each plot
  }
  \item{cex.lab}{\code{cex} of the labels
  }
  \item{verbose}{\code{Boolean} whether to write out extra information
  }
  \item{dev.printer}{Either \code{NULL} (default), which means that plotting is
    only to the screen, or a list with elements \code{device},
    \code{extension}, \code{width}, \code{height}, and \code{name}
    (with defaults \code{png}, \code{".png"}, \code{640}, \code{480},
    and \code{"Rplot"}, respectively) to \code{\link{dev.print}} the
    plots to files (only works in interactive sessions).
    This option is more or less obsolete since the \pkg{animation}
    package provides better features for output to files.
  }
  \item{\dots}{Extra arguments sent to the plot function.
  }
}

\author{
  Michael H\enc{ö}{oe}hle
}

\note{
  The \code{\link{animate.sts}} method provides a
  re-implementation and supersedes this function!
}

\seealso{
  Other \code{\link{stsplot}} types, and \code{\link{animate.sts}} for
  the new implementation.
}

\examples{
data("ha.sts")
print(ha.sts)

## map of total counts by district
plot(ha.sts, type=observed ~ 1 | unit)
## only show a sub-period total for two selected districts
plot(ha.sts[1:20,1:2], type=observed ~ 1 | unit)

\dontrun{
# space-time animation
plot(aggregate(ha.sts,nfreq=13), type= observed ~ 1 | unit * time)

#print the frames to a png device
#and do the animation without extra sleeping between frames
imgname <- file.path(tempdir(), "berlin")
plot(aggregate(ha.sts,nfreq=13), type = observed ~ 1 | unit * time,
     wait.ms=0, dev.printer=list(name=imgname))

#Use ImageMagick (you might have to adjust the path to 'convert')
system(paste0("convert -delay 50 ", imgname,
              "*.png ", imgname, "-animated.gif"))
}
}

\keyword{hplot}
\keyword{dynamic}
\keyword{spatial}
