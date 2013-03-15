\encoding{latin1}
\name{nowcast}
\alias{nowcast}
%Internal functions
\alias{nowcast.fit}
\alias{dist.median}
\alias{outside.ci}
\alias{logS}
\alias{RPS}

\title{
  Adjust observed epidemic curve for reporting delay of cases
}

\description{
  Nowcasting can help to obtain up-to-date information on epidemic
  trends during an outbreak situation, when reports about important
  indicators (such as hospitalizations) only arrive with delay.
  Note: This function is experimental and currently does not deal
  with right-truncation of the observations. An improved version is to
  be expected.
}

\usage{
nowcast(s, t, D, dEventCol = "dHospital", dReportCol = "dReport",
          method = c("freq.pi","bayes.nb","bayes.betapi","bayes.bnb","uniform"),
          aggregate.by="1 day",
          control = list(
            dRange = NULL,
            timeDelay = function(d1,d2) {as.numeric(d2-d1)},
            estimateF = NULL,
            alpha = 0.05,
            y.prior.max = 300,
            B = 1e5, 
            score = FALSE, PMF=FALSE, sts.truth=FALSE))
}

\arguments{
  \item{s}{
    an object of class \code{Date} denoting the day at which to do the
    nowcast
  }
  \item{t}{a vector of \code{Date} objects denoting the day(s) for which
    the projections are to be done. One needs to ensure t <= s.
  }
  \item{D}{A data frame with one row per case -- for each case on needs
  information on the day of the event (e.g. hospitalization) and the
  day of report of this event.
  }
  \item{dEventCol}{The name of the column in \code{D} which contains the
  date of the event, e.g. hospitalization. Default: \code{"dHospital"}.
  }
  \item{dReportCol}{Name of the column in \code{D} containing the date at
  which the report arrives at the respective register. Default:
  \code{"dReport"}.
  }
  \item{method}{A vector of strings denoting the different methods to
  use. Note that results of the first name in this list are officially
  returned by the function. However, it is possible to specify several
  methods here, e.g. in order to compare score evaluations.
  }
  \item{aggregate.by}{Time scale used for the temporal aggregation of
  the records in the data \code{D}. See 
  \link{linelist2sts} and \link{seq.Date} for further information.}
  \item{control}{A list with named arguments controlling the
    functionality of the nowcasting.
    \describe{
      \item{dRange}{Default: \code{NULL}. In this case the
	\code{dEventCol} column is used to extract the first and last
	available in \code{D}.}
      \item{timeDelay}{A function taking arguments \code{d1} and
    \code{d2} and which computes the time span between them. 
	Default: \code{function(d1,d2) {as.numeric(d2-d1)}}. However,
	more advanced function are imaginable. See examples.
      }
      \item{estimateF}{Denotes the estimation procedure for the CDF of
	the delay distribution. The function should have a header of the form
	\code{function(Ds,s,dReportCol)}, where \code{Ds} is the
        linelist of all cases available at \code{s} and \code{dReportCol}
        denotes the column containing the report data. If \code{NULL}
	(Default) then the empirical CDF is computed each time based on
        all available data. Otherwise, one can specify a function
	based on more complex modelling. Results needs to be a CDF, e.g.
	as a stepfunction.
      }
      \item{alpha}{Equal tailed (1-\eqn{\alpha}{alpha})*100\% prediction
	intervals are calculated. Default: 0.05.}
      \item{y.prior.max}{Limit of the support of \eqn{y_t}{yt}. The
	value needs to be high enough such that at this limit only little of
	the predictive distribution is right-truncated. Default: 300}
      \item{B=1e5}{Number of samples to use in the method
	\code{"freq.pi"}. Default: 1e5} 
      \item{score}{Compute scoring rules. Default: \code{FALSE}. The
  computed scores are found in the \code{control} slot of the result.}
      \item{PMF}{Boolean whether to teturn the probability mass
  functions of the individual forecasts (Default: \code{FALSE}). The
  result can be found in the \code{control} slot of the return object.}
\item{sts.truth}{Boolean whether to return an sts object containing the
  true number of counts (ignoring reporting delays). Result would be
  found in the \code{control} slot.}
      }
    }
  }
  
    

\details{
  No details currently available.    
}

\value{
  \code{nowcast} returns an object of \code{"\linkS4class{stsBP}"}. The
  \code{upperbound} slot contains the median of the method specified at
  position one of the argument \code{method}.
}

\references{
  Donker T, van Boven M, van Ballegooijen WM, Van't Klooster TM,
  Wielders CC, Wallinga J. (2011), Nowcasting pandemic influenza A/H1N1
  2009 hospitalizations in the Netherlands, Eur J Epidemiol. 2011
  Mar;26(3):195-201. Epub 2011 Mar 18.

  an der Heiden M, Wadl M, \enc{Höhle}{Hoehle} M (2011), Now-Casting
  during a huge outbreak of haemolytic-uremic syndrome in Germany 2011,
  Abstracts of the ESCAIDE conference, Stockholm, Sweden, 6-8 Nov, 2011.
}

\author{
Michael \enc{Höhle}{Hoehle}
}

\note{
  The methods and implementation is still experimental.
}

%\seealso{
%}

\examples{
#Load simulated outbreak data.
url <- paste("http://www.stat.uni-muenchen.de/~hoehle/",
             "teaching/moid2011/tutorials/cast-backnow/outbreak.txt",sep="")
D <- try(read.table(url,header=TRUE,colClasses=c("integer",rep("Date",3))))

if (!inherits(D, "try-error")) {
    s <- as.Date("2011-06-02") ;
    k <- 10 ; l <- 3
    t <- seq(s-k-l+1,s-l,by="1 day")
    dRange <- as.Date(c("2011-05-01","2011-07-10"))

    #Perform the nowcast
    nc <- nowcast(s=s,t=t,D=D,dEventCol="dHospital",dReportCol="dRegister",
                  method="bayes.nb",aggregate.by="1 day",control=list(dRange=dRange))

    #Show result
    plot(nc,xaxis.years=FALSE,dx.upperbound=0,legend=NULL,lty=c(1,1,1),
            lwd=c(1,1,2),ylab="Cases",xlab="Time (days)",main="",
            ylim=c(0,max(nc@ci,na.rm=TRUE)))
    idx <- max(which(!is.na(upperbound(nc))))
    lines( c(idx-0.5,idx+0.5), rep(upperbound(nc)[idx,],2),lwd=2,col="blue")

    ##Show CIs
    for (i in 1:nrow(nc)) {
      points(i, upperbound(nc)[i,], col="indianred")
      lines( i+c(-0.3,0.3), rep(nc@ci[i,,1],2),lty=1,col="indianred2")
      lines( i+c(-0.3,0.3), rep(nc@ci[i,,2],2),lty=1,col="indianred2")
      lines( rep(i,each=2), nc@ci[i,,],lty=2,col="indianred2")
    }
    #Add "now" on the x-axis
    points( as.numeric(s-dRange[1])+1,0,pch=10,cex=1.5,col="red")
}
}

\keyword{models}
\keyword{optimize}
