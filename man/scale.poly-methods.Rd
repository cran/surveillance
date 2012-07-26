\name{scale.poly,gpc.poly-method}
\docType{methods}
\alias{scale.poly,gpc.poly-method}
\title{Centering and Scaling a gpc.poly Polygon}
\description{
This is just a redefinition of the corresponding method from package
\pkg{gpclib} to also allow centering.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "gpc.poly")}}{
Shift (subtract) and scale (divide) the x and y coordinates of a \code{"gpc.poly"} object
by the amounts given in \code{center} and \code{scale}, respectively
(numeric vectors of length two).
}
}}
\keyword{methods}
\keyword{manip}

