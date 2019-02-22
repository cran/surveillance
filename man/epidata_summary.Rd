\name{epidata_summary}
\alias{summary.epidata}
\alias{print.summary.epidata}

\title{
Summarizing an Epidemic
}

\description{
  The \code{\link{summary}} method for \code{\link{class}}
  \code{"\link{epidata}"} gives an overview of the epidemic.  Its
  \code{\link{print}} method shows the type of the epidemic, the time range, the
  total number of individuals, the initially and never infected individuals and
  the size of the epidemic.  An excerpt of the returned \code{counters} data
  frame is also printed (see the Value section below).
}

\usage{
\method{summary}{epidata}(object, ...)

\method{print}{summary.epidata}(x, ...)
}

\arguments{
  \item{object}{an object inheriting from class \code{"epidata"}.}
  \item{x}{an object inheriting from class \code{"summary.epidata"}, i.e. an
    object returned by the function \code{summary.epidata}.}
  \item{\dots}{unused (argument of the generic).}
}

\value{
  A list with the following components:
  \item{type}{
    character string.  Compartmental type of the epidemic, i.e. one of "SIR",
    "SI", "SIS" or "SIRS".
  }
  \item{size}{
    integer.  Size of the epidemic, i.e. the number of initially susceptible
    individuals, which became infected during the course of the epidemic.
  }
  \item{initiallyInfected}{
    factor (with the same levels as the \code{id} column in the \code{"epidata"}
    object).  Set of initially infected individuals.
  }
  \item{neverInfected}{
    factor (with the same levels as the \code{id} column in the \code{"epidata"}
    object).  Set of never infected individuals, i.e. individuals, which were 
    neither initially infected nor infected during the course of the epidemic.
  }
  \item{coordinates}{
    numeric matrix of individual coordinates with as many rows as there are
    individuals and one column for each spatial dimension.  The row names of
    the matrix are the \code{id}s of the individuals.
  }
  \item{byID}{
    data frame with time points of infection and optionally removal and
    re-susceptibility (depending on the \code{type} of the epidemic) ordered
    by \code{id}.  If an event was not observed, the corresponding
    entry is missing.
  }
  \item{counters}{
    data frame containing all events (S, I and R) ordered by time.  The columns
    are \code{time}, \code{type} (of event), corresponding \code{id} and the
    three counters \code{nSusceptible}, \code{nInfectious} and \code{nRemoved}.
    The first row additionally shows the counters at the beginning of the
    epidemic, where the \code{type} and \code{id} column contain missing values.
  }
}

\author{
Sebastian Meyer
}

\seealso{
\code{\link{as.epidata}} for generating objects of class \code{"epidata"}.
}

\examples{
data("hagelloch")
s <- summary(hagelloch)
s          # uses the print method for summary.epidata
names(s)   # components of the list 's'

# positions of the individuals
plot(s$coordinates)

# events by id
head(s$byID)
}

\keyword{methods}
