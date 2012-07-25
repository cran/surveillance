\name{gpc.poly,coerce-methods}
\docType{methods}
\alias{coerce,gpc.poly,Polygons-method}
\alias{coerce,owin,gpc.poly-method}
\alias{coerce,Polygons,gpc.poly-method}
\alias{coerce,SpatialPolygons,gpc.poly-method}
\title{Some Additional Converters Between \pkg{gpclib}, \pkg{sp}, and \pkg{spatstat}}
\description{
The \pkg{surveillance} package defines some additional methods for conversion between
\code{"\link[gpclib:gpc.poly-class]{gpc.poly}"} (\pkg{gpclib}),
\code{"\link[sp:Polygons-class]{Polygons}"} (\pkg{sp}),
and \code{"\link[spatstat]{owin}"} (\pkg{spatstat}).
}
\section{Methods}{
\describe{
\item{\code{signature(from = "gpc.poly", to = "Polygons")}}{}
\item{\code{signature(from = "owin", to = "gpc.poly")}}{}
\item{\code{signature(from = "Polygons", to = "gpc.poly")}}{}
\item{\code{signature(from = "SpatialPolygons", to = "gpc.poly")}}{}
}}
\section{Details}{
Please note the varying polygon specifications in the packages:
\describe{
\item{\pkg{sp}:}{\emph{Repeat} first vertex at the end (closed),
anticlockwise order = hole, clockwise order = normal boundary}
\item{\pkg{spatstat}:}{do \emph{not repeat} first vertex,
anticlockwise order = normal boundary, clockwise order = hole}
\item{\pkg{gpclib}:}{\emph{not documented}, but it seems that it is prefered to
not repeat the first vertex and to have clockwise vertex order for the normal boundary
(and this is how we implemented the conversion here).}
}
}
\keyword{methods}
\keyword{spatial}

