%%% FIXME: remove this file and add all the aliases to sts-class.Rd ?
\name{stsSlot-generics}
\docType{methods}
\alias{alarms}
\alias{alarms<-}
\alias{upperbound}
\alias{upperbound<-}
\alias{control}
\alias{control<-}
\alias{epoch}
\alias{epoch<-}
\alias{observed}
\alias{observed<-}
\alias{population}
\alias{population<-}
\alias{multinomialTS}
\alias{multinomialTS<-}
\alias{neighbourhood}
\alias{neighbourhood<-}

\title{Generic functions to access \code{"sts"} slots}
\description{
 For almost every slot of the \code{"sts"} class, package \pkg{surveillance}
 defines a generic function of the same name (except for the \code{population}
 method where the slot is actually called \code{populationFrac}, and
 \code{alarms}, where the slot is actually called \code{alarm}) as well as
 a replacement version (\code{<-}) to extract or set the corresponding slot of
 a \code{sts} object.
 (This documentation is not really valid yet.)
}
%\section{Methods}{
%\describe{
%
%\item{x = "sts"}{ The slot of \code{x} is determined and returned or set.}
%}}
\seealso{
the \code{"\linkS4class{sts}"} class
}
\keyword{methods}
