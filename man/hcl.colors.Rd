\name{hcl.colors}
\alias{hcl.colors}
\title{
  HCL-based Heat Colors from the \pkg{colorspace} Package
}
\description{
  If package \pkg{colorspace} is available, its
  \link[colorspace]{heat_hcl} function is used to generate a color
  palette. Otherwise, the function simply calls
  \code{\link{heat.colors}} (or \code{\link{grey.colors}} if
  \code{use.color=FALSE}).
}
\usage{
hcl.colors(ncolors = 100, use.color = TRUE)
}
\arguments{
  \item{ncolors}{the number of colors (>= 1) to be in the palette.}
  \item{use.color}{logical. Should the palette use colors?
    Otherwise grey levels are returned.}
}
\value{
  A character vector of \code{ncolors} colors.
}
\seealso{
  \code{\link[colorspace]{heat_hcl}} in package \pkg{colorspace},
  or \code{\link{heat.colors}} in the base package \pkg{grDevices}.
}

\examples{
barplot(rep(1,10), col = hcl.colors(10), axes = FALSE)
}

\keyword{color}
\keyword{dplot}
