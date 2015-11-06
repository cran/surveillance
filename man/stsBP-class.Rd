\name{stsBP-class}
\Rdversion{1.1}
\docType{class}
\alias{stsBP-class}
% %
% % other names
% %\alias{aggregate,stsBP-method}
% %\alias{plot,stsBP-method}
% %\alias{plot,stsBP,missing-method}
% \alias{nrow,stsBP-method}
% \alias{ncol,stsBP-method}
% \alias{dim,stsBP-method}
% \alias{alarms,stsBP-method}
% \alias{alarms<-,stsBP-method}
% \alias{upperbound,stsBP-method}
% \alias{upperbound<-,stsBP-method}
% \alias{control,stsBP-method}
% \alias{control<-,stsBP-method}
% \alias{epoch,stsBP-method}
% \alias{epoch<-,stsBP-method}
% \alias{epochInYear,stsBP-method}
% \alias{observed,stsBP-method}
% \alias{observed<-,stsBP-method}
% \alias{population,stsBP-method}
% \alias{population<-,stsBP-method}
% \alias{multinomialTS,stsBP-method}
% \alias{multinomialTS<-,stsBP-method}
% \alias{neighbourhood,stsBP-method}
% \alias{neighbourhood<-,stsBP-method}
% \alias{colnames,stsBP,missing,missing-method}
% \alias{initialize,stsBP-method}
% \alias{as.data.frame,stsBP-method}

%Coerce method to convert to sts object
\alias{coerce,sts,stsBP-method}

\encoding{latin1}

\title{Class "stsBP" -- a class inheriting from class \code{sts} which
  allows the user to store the results of back-projecting or nowcasting
  surveillance time series} 

\description{
  A class inheriting from class \code{sts}, but with additional slots
  to store the result and associated confidence intervals from back
  projection of a \code{sts} object.

}

\section{Slots}{
  The slots are as for \code{"\linkS4class{sts}"}. However, two
  additional slots exists.
  
  \describe{
    \item{\code{ci}:}{An array containing the upper and lower limit
      of the confidence interval.}
    \item{\code{lambda}:}{Back projection component}
  }
}
\section{Methods}{
  The methods are the same as for \code{"\linkS4class{sts}"}.

  \itemize{
    \item{\code{signature(from = "sts", to = "stsBP")}}{
      Convert an object of class \code{sts} to class \code{stsBP}.
    }
  }
} %end section methods.


\author{M. \enc{Höhle}{Hoehle}}


\keyword{classes}
