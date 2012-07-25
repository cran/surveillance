\name{twinSIR_epidata}
\alias{as.epidata}
\alias{as.epidata.default}
\alias{as.epidata.epidataCS}
\alias{print.epidata}
\alias{[.epidata}
\alias{epidata}

\title{
Class for Epidemic Data Discrete in Space and Continuous in Time
}

\description{
  The function \code{as.epidata} is used to generate objects
  of \code{\link{class}} \code{"epidata"}.  Objects of this class are
  specific data frames containing the event history of an epidemic together
  with some additional attributes.  These objects are the basis for fitting
  spatio-temporal epidemic intensity models with the function
  \code{\link{twinSIR}}.  Note that the spatial information itself, i.e.
  the positions of the individuals, is assumed to be constant over time.  
  Besides epidemics following the SIR compartmental model, also data from SI,
  SIRS and SIS epidemics may be supplied.  Inference for the infectious process
  works as usual and simulation of such epidemics is also possible.
}

\usage{
as.epidata(data, ...)

\method{as.epidata}{default}(data, id.col, start.col, stop.col, atRiskY.col,
           event.col, Revent.col, coords.cols, f = list(), ...)
\method{as.epidata}{epidataCS}(data, tileCentroids, eps = 0.001, ...)

\method{print}{epidata}(x, ...)
\method{[}{epidata}(x, i, j, drop)
}

\arguments{
  \item{data}{
    For the default method a \code{\link{matrix}} or a \code{\link{data.frame}}.
    It contains the observed event history in a form similar to 
    \code{Surv(, type="counting")} with additional information (variables) along 
    the process.  It must not be sorted in any specific order; this will be done 
    automatically during conversion.  The observation period is splitted up into 
    \emph{consecutive} intervals of constant state - thus constant infection 
    intensities.  The data frame consists of a block of N (number of individuals) 
    rows for each of those time intervals (all rows in a block share the same start 
    and stop values... therefore the name \dQuote{block}), where there is one 
    row per individual in the block.  Each row describes the (fixed) state of 
    the individual during the interval given by the start and stop columns 
    \code{start.col} and \code{stop.col}.\cr
    Note that there may not be more than one event (infection or removal) in a
    single block.  Thus, in a single block, only one entry in the 
    \code{event.col} and \code{Revent.col} may be 1, all others are 0.  This
    rule follows the assumption that there are no concurrent events (infections
    or removals).
    
    Alternatively, \code{data} may be an object of class \code{"\link{epidataCS}"},
    i.e. space-continuous epidemic data used in the \code{\link{twinstim}} context,
    which will be converted (aggregated) to an object of class \code{"epidata"}
	for use with \code{\link{twinSIR}}.
  }
  \item{id.col}{
    single index of the \code{id} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).
    The \code{id} column identifies the individuals in the data frame.  It will
    be converted to a factor variable.
  }
  \item{start.col}{
    single index of the \code{start} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).
    The \code{start} column contains the (numeric) time points of the beginnings
    of the consecutive time intervals of the event history.  The minimum value
    in this column, i.e. the start of the observation period should be 0.
  }
  \item{stop.col}{
    single index of the \code{stop} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).
    The \code{stop} column contains the (numeric) time points of the ends
    of the consecutive time intervals of the event history.  The stop value must
    always be greater than the start value of a row.
  }
  \item{atRiskY.col}{
    single index of the \code{atRiskY} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).
    The \code{atRiskY} column indicates if the individual was \dQuote{at-risk}
    of becoming infected during the time interval (start; stop].  This variable 
    must be logical or in 0/1-coding.
    Individuals with \code{atRiskY == 0} in the first time interval (normally 
    the rows with \code{start == 0}) are taken as \emph{initially infectious}.
  }
  \item{event.col}{
    single index of the \code{event} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).
    The \code{event} column indicates if the individual became \emph{infected}
    at the \code{stop} time of the interval.  This variable must be logical or
    in 0/1-coding.
  }
  \item{Revent.col}{
    single index of the \code{Revent} column in \code{data}.  Can be numeric
    (by column number) or character (by column name).
    The \code{Revent} column indicates if the individual was \emph{recovered} 
    at the \code{stop} time of the interval.  This variable must be logical or
    in 0/1-coding.
  }
  \item{coords.cols}{
    index\emph{es} of the \code{coords} column\emph{s} in \code{data}. Can be
    a numeric (by column number) vector, a character (by column name) vector
    or \code{NULL} (in which case epidemic covariates are not calculateable).
    These columns contain the coordinates of the individuals.  It must be
    emphasized that the functions related to \code{\link{twinSIR}} models currently assume
    \emph{fixed positions} of the individuals during the whole epidemic.  Thus, an
    individual has the same coordinates in every block.  For simplicity, the
    coordinates are derived from the first time block only (normally the rows 
    with \code{start == 0}).
    The epidemic covariates are calculated based on the euclidian distances
    between the individuals, see \code{f}.
  }
  \item{f}{
    a \emph{named} list of distance functions or \code{list()} (the default),
    if calculation of epidemic covariates is not requested.  The functions must
    interact elementwise on a (distance) matrix so that - for a matrix D -
    f[[m]](D) results in a matrix.  A simple example is
    \code{function(u) {u <= 1}}, which indicates if the euclidian distance
    between the individuals is smaller than or equal to 1.  To ensure that an
    individual does not influence itself, the distance to itself is defined as
    \code{Inf}.  Consequently, all of the distance functions must have the
    property f[[m]](Inf) = 0.  The names of the functions will be the names of
    the epidemic variables in the resulting data frame.  The value of such a
    variable is computed as follows: \eqn{I(t)} denotes the set of infectious
    individuals just before time \eqn{t} and \eqn{s_i} the coordinate vector of
    individual \eqn{i}.  For individual \eqn{i} at time \eqn{t} the epidemic
    component \eqn{m} has the value
    \eqn{\sum_{j \in I(t)} f_m(||s_i - s_j||)}{%
         \sum_{j in I(t)} f[[m]](||s_i - s_j||)}
  }
  \item{tileCentroids}{
    a coordinate matrix of the tile centroids the row names of which are the tile levels.
    This will be the coordinates used for the population (i.e. the tiles from
    \code{"\link{epidataCS}"}) in the discrete-space \code{\link{twinSIR}} modelling.
  }
  \item{eps}{
  	numeric scalar for breaking tied removal and infection times between different
  	individuals (tiles), which might occur during conversion from
  	\code{"\link{epidataCS}"} to \code{"epidata"}. Rather dumb, this is simply done
  	by subtracting \code{eps} from each tied removal time.
  	One should consider other ways of breaking the tied event times.
  }
  \item{x}{
    an object of class \code{"epidata"}.
  }
  \item{\dots}{
    arguments passed to \code{\link{print.data.frame}}. Currently unused in \code{as.epidata}.
  }
  \item{i,j,drop}{
    arguments passed to \code{\link{[.data.frame}}.
  }
}

\details{
  The \code{print} method for objects of class \code{"epidata"} simply prints
  the data frame with a small header containing the time range of the observed
  epidemic and the number of infected individuals.  Usually, the data frames
  are quite long, so the summary method \code{\link{summary.epidata}} might be
  useful.  Also, indexing/subsetting \code{"epidata"} works exactly as for
  \code{\link[=[.data.frame]{data.frame}}s, but there is an own method, which
  assures consistency of the resulting \code{"epidata"} or drops this class, if
  necessary.
  
  SIS epidemics are implemented as SIRS epidemics where the length of the
  removal period equals 0.  This means that an individual, which has an R-event
  will be at risk immediately afterwards, i.e. in the following time block.
  Therefore, data of SIS epidemics have to be provided in that form containing
  \dQuote{pseudo-R-events}.
  
  Some comments on the conversion from \code{"\link{epidataCS}"} to \code{"epidata"}:
  the conversion results into SIS epidemics only,
  i.e. the at-risk indicator is set to 1 immediately after
  recovery. The lengths of the infectious periods are taken from
  \code{data$events$eps.t}. There will be no \code{f} columns in the resulting
  \code{"epidata"}. These must be generated by a subsequent call to
  \code{as.epidata} with desired \code{f}.
}

\note{
  The column name \code{"BLOCK"} is a reserved name.  This column will be 
  added automatically at conversion and the resulting data frame will be 
  sorted by this column and by id.  Also the names \code{"id"}, \code{"start"},
  \code{"stop"}, \code{"atRiskY"}, \code{"event"} and \code{"Revent"} are
  reserved for the respective columns only.
}

\value{
  a \code{data.frame} with the columns \code{"BLOCK"}, \code{"id"},
  \code{"start"}, \code{"stop"}, \code{"atRiskY"}, \code{"event"},
  \code{"Revent"} and the coordinate columns (with the original names from
  \code{data}), which are all obligatory.  These columns are followed by any 
  remaining columns of the input \code{data}.  Last but not least, the newly
  generated columns with epidemic variables corresponding to the functions
  in the list \code{f} are appended, if \code{length(f)} > 0.
  
  The \code{data.frame} is given the additional \emph{attributes}
  \item{"eventTimes"}{
    numeric vector of infection time points (sorted chronologically).
  }
  \item{"timeRange"}{
    numeric vector of length 2: \code{c(min(start), max(stop))}.
  }
  \item{"coords.cols"}{
    numeric vector containing the column indices of the coordinate columns in
    the resulting data frame.
  }
  \item{"f"}{
    this equals the argument \code{f}.
  }
}

\author{
Sebastian Meyer
}

\seealso{
The \code{\link[=plot.epidata]{plot}} and the
\code{\link[=summary.epidata]{summary}} method for class \code{"epidata"}.
Furthermore, the function \code{\link{animate.epidata}} for the animation of
epidemics.

Function \code{\link{twinSIR}} for fitting spatio-temporal epidemic intensity
models to epidemic data.

Function \code{\link{simEpidata}} for the simulation of epidemic data.
}

\examples{
# an artificial example of an event history from the package
data("foodata")
str(foodata)

# convert the data to an object of class "epidata",
# also generating some epidemic covariates
myEpidata <- as.epidata(foodata, id.col = 1, start.col = "start",
    stop.col = "stop", atRiskY.col = "atrisk", event.col = "infected",
    Revent.col = "removed", coords.cols = c("x","y"),
    f = list(B1 = function(u) u<=1,
             B2 = function(u) u>1 & is.finite(u))
)
# note the is.finite restriction in B2 to ensure that f[[i]](Inf) = 0, for all i

str(myEpidata)
subset(myEpidata, BLOCK == 1)

summary(myEpidata)          # see 'summary.epidata'
plot(myEpidata)             # see 'plot.epidata' and also 'animate.epidata'
stateplot(myEpidata, "15")  # see 'stateplot'

\dontrun{ # works in interactive mode, but not in R CMD check
data("fooepidata")
stopifnot(identical(myEpidata, fooepidata))
}
}

\keyword{spatial}
\keyword{classes}
\keyword{manip}
