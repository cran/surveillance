\name{epidata_animate}
\alias{animate.epidata}
\alias{animate.summary.epidata}

\title{
Spatio-Temporal Animation of an Epidemic
}

\description{
  Function for the animation of epidemic data, i.e. objects inheriting from 
  class \code{"epidata"}.  This only works with 1- or 2-dimensional coordinates 
  and is not useful if some individuals share the same coordinates 
  (overlapping).  There are two types of animation, see argument 
  \code{time.spacing}.  Besides the direct plotting in the \R session, it is 
  also possible to generate a sequence of graphics files to create animations 
  outside \R.
}

\usage{
\method{animate}{summary.epidata}(object, main = "An animation of the epidemic",
        pch = 19, col = c(3, 2, gray(0.6)), time.spacing = NULL,
        sleep = quote(5/.nTimes), legend.opts = list(), timer.opts = list(),
        end = NULL, generate.snapshots = NULL, ...)

\method{animate}{epidata}(object, ...)
}

\arguments{
  \item{object}{
    an object inheriting from class \code{"epidata"} or
    \code{"summary.epidata"}.  In the former case, its summary is calculated
    and the function continues as in the latter case, passing all \code{...}
    arguments to the \code{summary.epidata} method.
  }
  \item{main}{
    a main title for the plot, see also \code{\link{title}}.
  }
  \item{pch, col}{
    vectors of length 3 specifying the point symbols and colors for
    susceptible, infectious and removed individuals (in this order).
    The vectors are recycled if necessary.
    By default, susceptible individuals are marked as filled green circles,
    infectious individuals as filled red circles and removed individuals as
    filled gray circles.  Note that the symbols are iteratively drawn
    (overlayed) in the same plotting region as time proceeds.
    For information about the possible values of \code{pch} and \code{col}, see
    the help pages of \code{\link{points}} and \code{\link{par}}, respectively.
  }
  \item{time.spacing}{
    time interval for the animation steps.  If \code{NULL} (the default), the
    events are plotted one by one with pauses of \code{sleep} seconds.  Thus, 
    it is just the \emph{ordering} of the events, which is shown.  To plot
    the appearance of events proportionally to the exact time line,
    \code{time.spacing} can be set to a numeric value indicating the period of
    time between consecutive plots.  Then, for each time point in
    \code{seq(0, end, by = time.spacing)} the current state of the epidemic can
    be seen and an additional timer indicates the current time
    (see \code{timer.opts} below).  The argument \code{sleep} will be the
    artificial pause in seconds between two of those time points.
  }
  \item{sleep}{
    time in seconds to \code{\link{Sys.sleep}} before the next plotting
    event.  By default, each artificial pause is of length \code{5/.nTimes}
    seconds, where \code{.nTimes} is the number of events (infections and
    removals) of the epidemic, which is evaluated in the function body.
    Thus, for \code{time.spacing = NULL} the animation has a duration of
    approximately 5 seconds. In the other case, \code{sleep} is the duration of
    the artificial pause between two time points.
    Note that \code{sleep} is ignored on non-interactive devices
    (see \code{\link{dev.interactive}})
  }
  \item{legend.opts}{
    either a list of arguments passed to the \code{\link{legend}} function or
    \code{NULL} (or \code{NA}), in which case no legend will be plotted.  All 
    necessary arguments have sensible defaults and need not be specified, i.e.
    \describe{
      \item{\code{x}:}{\code{"topright"}}
      \item{\code{legend}:}{\code{c("susceptible", "infectious", "removed")}}
      \item{\code{pch}:}{same as argument \code{pch} of the main function}
      \item{\code{col}:}{same as argument \code{col} of the main function}
    }
  }
  \item{timer.opts}{
    either a list of arguments passed to the \code{\link{legend}} function or
    \code{NULL} (or \code{NA}), in which case no timer will be plotted.  All 
    necessary arguments have sensible defaults and need not be specified, i.e.
    \describe{
      \item{\code{x}:}{\code{"bottomright"}}
      \item{\code{title}:}{\code{"time"}}
      \item{\code{box.lty}:}{\code{0}}
      \item{\code{adj}:}{\code{c(0.5,0.5)}}
      \item{\code{inset}:}{\code{0.01}}
      \item{\code{bg}:}{\code{"white"}}
    }
    Note that the argument \code{legend}, which is the current time of the
    animation, can not be modified.
  }
  \item{end}{
    ending time of the animation in case of \code{time.spacing} not being
    \code{NULL}.  By default (\code{NULL}), time stops after the last event.
  }
  \item{generate.snapshots}{
    By default (\code{NULL}), the animation is not saved to image files
    but only shown on the on-screen device. In order to print to files,
    \code{time.spacing} must not be \code{NULL}, a screen device must be
    available, and there are two options:\cr
    If the framework of the \pkg{animation} package
    should be used, i.e. the \code{animate}-call is passed as the
    \code{expr} argument to one of the \code{save*} functions of the
    \pkg{animation} package, then set
    \code{generate.snapshots = img.name}, where \code{img.name} is the
    base name for the generated images (the same as passed to the
    \code{save*} function). The path and format (type, width, height)
    for the generated images is derived from
    \code{\link[animation]{ani.options}}. See the last example below.\cr
    Alternatively, \code{generate.snapshots} may be a list of arguments
    passed to the function \code{\link{dev.print}}, which then is
    executed at each time point of the grid defined by
    \code{time.spacing}.  Essentially, this is used for 
    saving the produced snapshots to files, e.g.

    \code{generate.snapshots = %
    list(device=pdf, file=quote(paste("epidemic_",sprintf(form,tp),".pdf",%
    sep="")))}

    will store the animation steps in pdf-files in the current
    working directory, where the file names each end with the time point
    represented by the corresponding plot.  Because the variables \code{tp}
    and \code{form} should only be evaluated inside the function the
    \code{file} argument is \code{quote}d. Alternatively, the file name
    could also make use of the internal plot index \code{i}, e.g., use
    \code{file=quote(paste("epidemic",i,".pdf",sep=""))}.
  }
  \item{\dots}{
    further graphical parameters passed to the basic call of \code{plot}, e.g.
    \code{las}, \code{cex.axis} (etc.) and \code{mgp}.
  }
}

%\value{
%  invisibly returns \code{NULL}.
%}

\author{
Sebastian Meyer
}

\seealso{
\code{\link{summary.epidata}} for the data, on which the plot is based.
\code{\link{plot.epidata}} for plotting the evolution of an epidemic by
the numbers of susceptible, infectious and removed individuals.

The contributed \R package \pkg{animation}.
}

\examples{
data("hagelloch")
(s <- summary(hagelloch))

# plot the ordering of the events only
animate(s)   # or: animate(hagelloch)

# with timer (animate only up to t=10)
animate(s, time.spacing=0.1, end=10, sleep=0.01,
        legend.opts=list(x="topleft"))

# Such an animation can be saved in various ways using tools of
# the animation package, e.g., saveHTML()
if (interactive() && require("animation")) {
  oldwd <- setwd(tempdir())  # to not clutter up the current working dir
  saveHTML({
    par(bg="white")  # default "transparent" is grey in some browsers
    animate(s, time.spacing=1, sleep=0, legend.opts=list(x="topleft"),
            generate.snapshots="epiani")
  }, use.dev=FALSE, img.name="epiani", ani.width=600, interval=0.5)
  setwd(oldwd)
}
}

\keyword{hplot}
\keyword{dynamic}
\keyword{spatial}
