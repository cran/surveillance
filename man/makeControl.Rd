\name{makeControl}
\alias{makeControl}
\title{Generate \code{control} Settings for an \code{hhh4} Model}
\usage{
makeControl(f = list(~1), S = list(0, 0, 1), period = 52, offset = 1, ...)
}
\arguments{
\item{f, S, period}{
  arguments for \code{\link{addSeason2formula}} defining
  each of the three model formulae in the order (\code{ar}, \code{ne},
  \code{end}). Recycled if necessary within \code{\link{mapply}}.
}
\item{offset}{
  multiplicative component offsets in the order (\code{ar},
  \code{ne}, \code{end}).
}
\item{...}{
  further elements for the \code{\link{hhh4}} control list. The
  \code{family} parameter is set to \code{"NegBin1"} by default.
}
}
\value{
a list for use as the \code{control} argument in \code{\link{hhh4}}.
}
\description{
Generate \code{control} Settings for an \code{hhh4} Model
}
\examples{
makeControl()

## a simplistic model for the fluBYBW data
## (first-order transmission only, no district-specific intercepts)
data("fluBYBW")
mycontrol <- makeControl(
    f = list(~1, ~1, ~t), S = c(1, 1, 3),
    offset = list(population(fluBYBW)),  # recycled -> in all components
    ne = list(normalize = TRUE), verbose = TRUE)
str(mycontrol)
\dontrun{fit <- hhh4(fluBYBW, mycontrol)}
}
\author{
Sebastian Meyer
}
