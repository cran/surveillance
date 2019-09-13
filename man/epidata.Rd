\encoding{latin1}
\name{epidata}
\alias{as.epidata}
\alias{as.epidata.data.frame}
\alias{as.epidata.default}
\alias{print.epidata}
\alias{[.epidata}
\alias{update.epidata}
\alias{epidata}

\title{
Continuous-Time SIR Event History of a Fixed Population
}

\description{
  The function \code{as.epidata} is used to generate objects
  of class \code{"epidata"}.  Objects of this class are
  specific data frames containing the event history of an epidemic together
  with some additional attributes.  These objects are the basis for fitting
  spatio-temporal epidemic intensity models with the function
  \code{\link{twinSIR}}.  Their implementation is illustrated in Meyer
  et al. (2017, Section 4), see \code{vignette("twinSIR")}.
  Note that the spatial information itself, i.e.
  the positions of the individuals, is assumed to be constant over time.  
  Besides epidemics following the SIR compartmental model, also data from SI,
  SIRS and SIS epidemics may be supplied.
}

\usage{
as.epidata(data, ...)

\method{as.epidata}{data.frame}(data, t0,
           tE.col, tI.col, tR.col, id.col, coords.cols,
           f = list(), w = list(), D = dist,
           max.time = NULL, keep.cols = TRUE, ...)
\method{as.epidata}{default}(data, id.col, start.col, stop.col,
           atRiskY.col, event.col, Revent.col, coords.cols,
           f = list(), w = list(), D = dist, .latent = FALSE, ...)

\method{print}{epidata}(x, ...)
\method{[}{epidata}(x, i, j, drop)
\method{update}{epidata}(object, f = list(), w = list(), D = dist, ...)
}

\arguments{
  \item{data}{
    For the \code{data.frame}-method, a data frame with as many rows as
    there are individuals in the population and time columns indicating
    when each individual became exposed (optional), infectious
    (mandatory, but can be \code{NA} for non-affected individuals) and
    removed (optional). Note that this data format does not allow for
    re-infection (SIRS) and time-varying covariates.
    The \code{data.frame}-method converts the individual-indexed data
    frame to the long event history start/stop format and then feeds it
    into the default method. If calling the generic function
    \code{as.epidata} on a \code{data.frame} and the \code{t0} argument
    is missing, the default method is called directly.\cr
    For the default method, \code{data} can be a \code{\link{matrix}} or
    a \code{\link{data.frame}}.
    It must contain the observed event history in a form similar to 
    \code{Surv(, type="counting")} in package \pkg{survival},
    with additional information (variables) along 
    the process.  Rows will be sorted automatically during conversion.
    The observation period is split up into \emph{consecutive}
    intervals of constant state - thus constant infection intensities.
    The data frame consists of a block of \eqn{N} (number of individuals) 
    rows for each of those time intervals (all rows in a block have the same start 
    and stop values\dots therefore the name \dQuote{block}), where there is one 
    row per individual in the block.  Each row describes the (fixed) state of 
    the individual during the interval given by the start and stop columns 
    \code{start.col} and \code{stop.col}.\cr
    Note that there may not be more than one event (infection or removal) in a
    single block.  Thus, in a single block, only one entry in the 
    \code{event.col} and \code{Revent.col} may be 1, all others are 0.  This
    rule follows the point process characteristic that there are no
    concurrent events (infections or removals).
  }
  \item{t0,max.time}{
    observation period. In the resulting \code{"epidata"}, the time
    scale will be relative to the start time \code{t0}.
    Individuals that have already been removed prior to \code{t0}, i.e.,
    rows with \code{tR <= t0}, will be dropped.
    The end of the observation period (\code{max.time}) will by default
    (\code{NULL}, or if \code{NA}) coincide with the last observed event.
  }
  \item{tE.col, tI.col, tR.col}{
    single numeric or character indexes of the time columns in
    \code{data}, which specify when the individuals became exposed,
    infectious and removed, respectively.
    \code{tE.col} and \code{tR.col} can be missing, corresponding to
    SIR, SEI, or SI data. \code{NA} entries mean that the respective
    event has not (yet) occurred. Note that \code{is.na(tE)} implies
    \code{is.na(tI)} and \code{is.na(tR)}, and \code{is.na(tI)} implies
    \code{is.na(tR)} (and this is checked for the provided data).\cr
    CAVE: Support for latent periods (\code{tE.col}) is experimental!
          \code{\link{twinSIR}} cannot handle them anyway.
  }
  \item{id.col}{
    single numeric or character index of the \code{id} column in \code{data}.
    The \code{id} column identifies the individuals in the data frame.
    It is converted to a factor by calling \code{\link{factor}}, i.e.,
    unused levels are dropped if it already was a factor.
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
    numeric (by column number), character (by column name), or \code{NULL}
    (no coordinates, e.g., if \code{D} is a pre-specified distance matrix).
    These columns contain the individuals' coordinates, which determine
    the distance matrix for the distance-based components of the force
    of infection (see argument \code{f}). By default, Euclidean distance
    is used (see argument \code{D}).\cr
    Note that the functions related to \code{\link{twinSIR}} currently assume
    \emph{fixed positions} of the individuals during the whole epidemic.  Thus, an
    individual has the same coordinates in every block.  For simplicity, the
    coordinates are derived from the first time block only (normally the rows 
    with \code{start == 0}).\cr
    The \code{\link[=animate.epidata]{animate}}-method requires coordinates.
  }
  \item{f}{
    a \emph{named} list of \emph{vectorized} functions for a
    distance-based force of infection.
    The functions must interact elementwise on a (distance) matrix \code{D} so that
    \code{f[[m]](D)} results in a matrix.  A simple example is
    \code{function(u) {u <= 1}}, which indicates if the Euclidean distance
    between the individuals is smaller than or equal to 1.
    The names of the functions determine the names of the epidemic variables
    in the resulting data frame.  So, the names should not coincide with
    names of other covariates.
    The distance-based weights are computed as follows:
    Let \eqn{I(t)} denote the set of infectious
    individuals just before time \eqn{t}.
    Then, for individual \eqn{i} at time \eqn{t}, the
    \eqn{m}'th covariate has the value
    \eqn{\sum_{j \in I(t)} f_m(d_{ij})}{%
         \sum_{j in I(t)} f[[m]](d[i,j])},
    where \eqn{d_{ij}}{d[i,j]} denotes entries of the distance matrix
    (by default this is the Euclidean distance \eqn{||s_i - s_j||}
    between the individuals' coordinates, but see argument \code{D}).
  }
  \item{w}{
    a \emph{named} list of \emph{vectorized} functions for extra 
    covariate-based weights \eqn{w_{ij}}{w_ij} in the epidemic component.
    Each function operates on a single time-constant covariate in
    \code{data}, which is determined by the name of the first argument:
    The two function arguments should be named \code{varname.i} and
    \code{varname.j}, where \code{varname} is one of \code{names(data)}.
    Similar to the components in \code{f}, \code{length(w)} epidemic
    covariates will be generated in the resulting \code{"epidata"} named
    according to \code{names(w)}.  So, the names should not coincide with
    names of other covariates.  For individual \eqn{i} at time
    \eqn{t}, the \eqn{m}'th such covariate has the value
    \eqn{\sum_{j \in I(t)} w_m(z^{(m)}_i, z^{(m)}_j)},
    where \eqn{z^{(m)}} denotes the variable in \code{data} associated
    with \code{w[[m]]}.
  }
  \item{D}{
    either a function to calculate the distances between the individuals
    with locations taken from \code{coord.cols} (the default is
    Euclidean distance via the function \code{\link{dist}}) and
    the result converted to a matrix via \code{\link{as.matrix}},
    or a pre-computed distance matrix with \code{dimnames} containing
    the individual ids (a classed \code{"\linkS4class{Matrix}"} is supported).
  }
  \item{keep.cols}{
    logical indicating if all columns in \code{data}
    should be retained (and not only the obligatory \code{"epidata"}
    columns), in particular any additional columns with 
    time-constant individual-specific covariates.
    Alternatively, \code{keep.cols} can be a numeric or character vector
    indexing columns of \code{data} to keep.
  }
  \item{.latent}{
    (internal) logical indicating whether to allow for latent periods
    (EXPERIMENTAL). Otherwise (default), the function verifies that an
    event (i.e., switching to the I state) only happens when the
    respective individual is at risk (i.e., in the S state).
  }
  \item{x,object}{
    an object of class \code{"epidata"}.
  }
  \item{\dots}{
    arguments passed to \code{\link{print.data.frame}}. Currently unused
    in the \code{as.epidata}-methods.
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
  The \code{update}-method can be used to add or replace distance-based
  (\code{f}) or covariate-based (\code{w}) epidemic variables in an
  existing \code{"epidata"} object.
  
  SIS epidemics are implemented as SIRS epidemics where the length of the
  removal period equals 0.  This means that an individual, which has an R-event
  will be at risk immediately afterwards, i.e. in the following time block.
  Therefore, data of SIS epidemics have to be provided in that form containing
  \dQuote{pseudo-R-events}.  
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
  \item{"w"}{
    this equals the argument \code{w}.
  }
}

\author{
Sebastian Meyer
}

\seealso{
The \code{\link{hagelloch}} data as an example.

The \code{\link[=plot.epidata]{plot}} and the
\code{\link[=summary.epidata]{summary}} method for class \code{"epidata"}.
Furthermore, the function \code{\link{animate.epidata}} for the animation of
epidemics.

Function \code{\link{twinSIR}} for fitting spatio-temporal epidemic intensity
models to epidemic data.

Function \code{\link{simEpidata}} for the simulation of epidemic data.
}

\references{
  Meyer, S., Held, L. and \enc{Höhle}{Hoehle}, M. (2017):
  Spatio-temporal analysis of epidemic phenomena using the \R package
  \pkg{surveillance}.
  \emph{Journal of Statistical Software}, \bold{77} (11), 1-55.
  \doi{10.18637/jss.v077.i11}
}

\examples{
data("hagelloch")   # see help("hagelloch") for a description
head(hagelloch.df)

## convert the original data frame to an "epidata" event history
myEpi <- as.epidata(hagelloch.df, t0 = 0,
                    tI.col = "tI", tR.col = "tR", id.col = "PN",
                    coords.cols = c("x.loc", "y.loc"),
                    keep.cols = c("SEX", "AGE", "CL"))
\dontshow{
## test consistency with default method
evHist <- as.data.frame(myEpi)[,-1]
myEpi2 <- as.epidata(
    evHist, id.col = 1, start.col = "start", stop.col = "stop",
    atRiskY.col = "atRiskY", event.col = "event", Revent.col = "Revent",
    coords.cols = c("x.loc", "y.loc")
)
stopifnot(identical(myEpi, myEpi2))
}

str(myEpi)
head(as.data.frame(myEpi))  # "epidata" has event history format
summary(myEpi)              # see 'summary.epidata'
plot(myEpi)                 # see 'plot.epidata' and also 'animate.epidata'

## add distance- and covariate-based weights for the force of infection
## in a twinSIR model, see vignette("twinSIR") for a description
myEpi <- update(myEpi,
    f = list(
        household    = function(u) u == 0,
        nothousehold = function(u) u > 0
    ),
    w = list(
        c1 = function (CL.i, CL.j) CL.i == "1st class" & CL.j == CL.i,
        c2 = function (CL.i, CL.j) CL.i == "2nd class" & CL.j == CL.i
    )
)
## this is now identical to the prepared hagelloch "epidata"
stopifnot(all.equal(myEpi, hagelloch))

\dontshow{
## test with precomputed distance matrix D
myEpi3 <- suppressWarnings( # from overwriting existing f columns
    update(hagelloch, f = attr(hagelloch, "f"),
           D = as.matrix(dist(hagelloch.df[c("x.loc", "y.loc")])))
)
stopifnot(identical(hagelloch, myEpi3))
}
}

\keyword{spatial}
\keyword{classes}
\keyword{manip}
