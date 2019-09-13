\name{hcl.colors}
\alias{.hcl.colors}
\title{
  HCL-based Heat Colors from the \pkg{colorspace} Package
}
\description{
  If package \pkg{colorspace} is available, its
  \link[colorspace]{heat_hcl} function is used to generate a color
  palette. Otherwise, the similar Heat 2 palette from R's own
  \code{hcl.colors} (in R >= 3.6.0) or basic
  \code{\link{heat.colors}} are used.

  This function was exported as \code{hcl.colors} in \pkg{surveillance}
  1.14.0 - 1.17.0 but is now internal to avoid a name clash with R 3.6.0
  (or later), which introduced a function of that name in the base
  package \pkg{grDevices}.
}
\usage{
.hcl.colors(ncolors = 100, use.color = TRUE)
}
\arguments{
  \item{ncolors}{the number of colors (>= 1) to be in the palette.}
  \item{use.color}{logical. Should the palette use colors?
    Otherwise grey levels are returned.}
}
\value{
  A character vector of \code{ncolors} colors.
}

\examples{
barplot(rep(1,10), col = surveillance:::.hcl.colors(10), axes = FALSE)
}

\keyword{color}
\keyword{dplot}
\keyword{internal}
