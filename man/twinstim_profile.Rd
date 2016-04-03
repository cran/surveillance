\encoding{latin1}
\name{twinstim_profile}
\alias{profile.twinstim}

\title{
Profile Likelihood Computation and Confidence Intervals for
\code{twinstim} objects
}

\description{
  Function to compute estimated and profile likelihood based confidence
  intervals for \code{twinstim} objects. Computations might be cumbersome!
  
  WARNING: the implementation is not well tested, simply uses
  \code{optim} (ignoring optimizer settings from the original fit),
  and does not return the complete set of coefficients at each grid point.
}

\usage{
\method{profile}{twinstim}(fitted, profile, alpha = 0.05,
        control = list(fnscale = -1, maxit = 100, trace = 1),
        do.ltildeprofile=FALSE, ...)
}

\arguments{
  \item{fitted}{
    an object of class \code{"twinstim"}.
  }
  \item{profile}{
    a list with elements being numeric vectors of length 4.  These vectors must
    have the form \code{c(index, lower, upper, gridsize)}.
    \describe{
      \item{\code{index}:}{
        index of the parameter to be profiled in the vector \code{coef(fitted)}.
      }
      \item{\code{lower, upper}:}{
        lower/upper limit of the grid on which the profile log-likelihood is
        evaluated. Can also be \code{NA} in which case \code{lower/upper} equals 
        the lower/upper bound of the respective 0.3 \% Wald confidence interval
        (+-3*se).
      }
      \item{\code{gridsize}:}{
        grid size of the equally spaced grid between lower and upper.  Can also
        be 0 in which case the profile log-likelihood for this parameter is not 
        evaluated on a grid.
      }
    }
  }
  \item{alpha}{
    \eqn{(1-\alpha)\%}{(1-alpha)\%} profile likelihood based confidence
    intervals are computed.  If alpha <= 0, then no confidence intervals are
    computed. This is currently not implemented.
  }
  \item{control}{
    control object to use in \code{\link{optim}} for the profile log-likelihood
    computations. It might be necessary to control \code{maxit} or
  \code{reltol} in order to obtain results in finite time.
  }
  \item{do.ltildeprofile}{If \code{TRUE} calculate profile likelihood as
    well. This might take a while, since an optimisation for all other
    parameters has to be performed. Useful for likelihood based
    confidence intervals. Default: \code{FALSE}.
  }
  \item{\dots}{
    unused (argument of the generic).
  }
}

\value{
list with profile log-likelihood evaluations on the grid, and
-- not implemented yet --
highest likelihood and Wald confidence intervals.
The argument \code{profile} is also returned.
}

\author{
Michael \enc{Höhle}{Hoehle}
}

\examples{
# profiling takes a while
\dontrun{
#Load the twinstim model fitted to the IMD data
data("imdepi", "imdepifit")
# for profiling we need the model environment
imdepifit <- update(imdepifit, model=TRUE)

#Generate profiling object for a list of parameters for the new model
names <- c("h.(Intercept)","e.typeC")
coefList <- lapply(names, function(name) {
  c(pmatch(name,names(coef(imdepifit))),NA,NA,11)
})

#Profile object (necessary to specify a more loose convergence
#criterion). Speed things up by using do.ltildeprofile=FALSE (the default)
prof <- profile(imdepifit, coefList,
  control=list(reltol=0.1, REPORT=1), do.ltildeprofile=TRUE)

#Plot result for one variable
par(mfrow=c(1,2))
for (name in names) {
  with(as.data.frame(prof$lp[[name]]),
       matplot(grid,cbind(profile,estimated,wald),
               type="l",xlab=name,ylab="loglik"))
  legend(x="bottomleft",c("profile","estimated","wald"),lty=1:3,col=1:3)
}
}
}

\keyword{htest}
\keyword{methods}
\keyword{optimize}
\keyword{dplot}
