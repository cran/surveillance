\name{pairedbinCUSUM}
\alias{pairedbinCUSUM}
\alias{pairedbinCUSUM.runlength}
\alias{pairedbinCUSUM.LLRcompute}
\encoding{latin1}

\title{Paired binary CUSUM and its run-length computation}
\description{
  CUSUM for paired binary data as described in Steiner et al. (1999).
}
\usage{
pairedbinCUSUM(stsObj, control = list(range=NULL,theta0,theta1,
                                      h1,h2,h11,h22))
pairedbinCUSUM.runlength(p,w1,w2,h1,h2,h11,h22, sparse=FALSE)
}
\arguments{
  \item{stsObj}{Object of class \code{sts} containing the paired
    responses for each of the, say n, patients. The observed slot of
    \code{stsObj} is thus a \eqn{n \times 2}{n x 2} matrix.}
  \item{control}{Control object as a list containing several parameters.
    \itemize{
      \item{\code{range}}{Vector of indices in the observed slot to monitor.}
      \item{\code{theta0}}{In-control parameters of the paired binary CUSUM.}
      \item{\code{theta1}}{Out-of-control parameters of the paired binary CUSUM.}
      \item{\code{h1}}{Primary control limit (=threshold) of 1st CUSUM.}
      \item{\code{h2}}{Primary control limit (=threshold) of 2nd CUSUM.}
      \item{\code{h11}}{Secondary limit for 1st CUSUM.}
      \item{\code{h22}}{Secondary limit for 2nd CUSUM.}
    }
  }
  \item{p}{Vector giving the probability of the four different
    possible states, i.e. c((death=0,near-miss=0),(death=1,near-miss=0),
       (death=0,near-miss=1),(death=1,near-miss=1)).}
  \item{w1}{The parameters \code{w1} and \code{w2} are the sample
    weights vectors for the two CUSUMs, see eqn. (2) in the paper. We
    have that \code{w1} is equal to deaths }
  \item{w2}{As for \code{w1}}
  \item{h1}{decision barrier for 1st individual cusums}
  \item{h2}{decision barrier for 2nd cusums}
  \item{h11}{together with \code{h22} this makes up the joing decision barriers}
  \item{h22}{together with \code{h11} this makes up the joing decision barriers}
  \item{sparse}{Boolean indicating whether to use sparse matrix
    computations from the \code{Matrix} library (usually much faster!). Default:
  \code{FALSE}.}
}
\details{
  For details about the method see the Steiner et al. (1999) reference
  listed below. Basically, two individual CUSUMs are run each based on a
  logistic regression model. The combined CUSUM not only signals if one
  of its two individual CUSUMs signals, but also if the two CUSUMs
  simultaneously cross the secondary limits.
}
\seealso{\code{\link{categoricalCUSUM}}}
\value{An \code{sts} object with \code{observed}, \code{alarm},
  etc. slots trimmed to the \code{control$range} indices.
}
\references{
  Steiner, S. H., Cook, R. J., and Farewell, V. T. (1999), Monitoring
  paired binary surgical outcomes using cumulative sum charts,
  Statistics in Medicine, 18, pp. 69--86.
}
\examples{
#Set in-control and out-of-control parameters as in paper
theta0 <- c(-2.3,-4.5,2.5)
theta1 <- c(-1.7,-2.9,2.5)

#Small helper function to compute the paired-binary likelihood
#of the length two vector yz when the true parameters are theta
dPBin <- function(yz,theta) {
    exp(dbinom(yz[1],size=1,prob=plogis(theta[1]),log=TRUE) +
    dbinom(yz[2],size=1,prob=plogis(theta[2]+theta[3]*yz[1]),log=TRUE))
}

#Likelihood ratio for all four possible configurations
p <- c(dPBin(c(0,0), theta=theta0), dPBin(c(0,1), theta=theta0),
       dPBin(c(1,0), theta=theta0), dPBin(c(1,1), theta=theta0))

#Compute ARL using non-sparse matrix operations
\dontrun{
pairedbinCUSUM.runlength(p,w1=c(-1,37,-9,29),w2=c(-1,7),h1=70,h2=32,h11=38,h22=17)
}

#Sparse computations don't work on all machines (e.g. the next line
#might lead to an error. If it works this call can be considerably (!) faster
#than the non-sparse call.
\dontrun{
pairedbinCUSUM.runlength(p,w1=c(-1,37,-9,29),w2=c(-1,7),h1=70,h2=32,
                         h11=38,h22=17,sparse=TRUE)
}

#Use paired binary CUSUM on the De Leval et al. (1994) arterial switch
#operation data on 104 newborn babies
data("deleval")

#Switch between death and near misses
observed(deleval) <- observed(deleval)[,c(2,1)]

#Run paired-binary CUSUM without generating alarms.
pb.surv <- pairedbinCUSUM(deleval,control=list(theta0=theta0,
             theta1=theta1,h1=Inf,h2=Inf,h11=Inf,h22=Inf))

plot(pb.surv, xaxis.labelFormat=NULL, ylab="CUSUM Statistic")



######################################################################
#Scale the plots so they become comparable to the plots in Steiner et
#al. (1999). To this end a small helper function is defined.
######################################################################

######################################################################
#Log LR for conditional specification of the paired model
######################################################################
LLR.pairedbin <- function(yz,theta0, theta1) {
    #In control
    alphay0 <- theta0[1] ; alphaz0 <- theta0[2] ; beta0 <- theta0[3]
    #Out of control
    alphay1 <- theta1[1] ; alphaz1 <- theta1[2] ; beta1 <- theta1[3]
    #Likelihood ratios
    llry <- (alphay1-alphay0)*yz[1]+log(1+exp(alphay0))-log(1+exp(alphay1))
    llrz <- (alphaz1-alphaz0)*yz[2]+log(1+exp(alphaz0+beta0*yz[1]))-
                                    log(1+exp(alphaz1+beta1*yz[1]))
    return(c(llry=llry,llrz=llrz))
}


val <- expand.grid(0:1,0:1)
table <- t(apply(val,1, LLR.pairedbin, theta0=theta0, theta1=theta1))
w1 <- min(abs(table[,1]))
w2 <- min(abs(table[,2]))
S <- upperbound(pb.surv) / cbind(rep(w1,nrow(observed(pb.surv))),w2)

#Show results
par(mfcol=c(2,1))
plot(1:nrow(deleval),S[,1],type="l",main="Near Miss",xlab="Patient No.",
     ylab="CUSUM Statistic")
lines(c(0,1e99), c(32,32),lty=2,col=2)
lines(c(0,1e99), c(17,17),lty=2,col=3)

plot(1:nrow(deleval),S[,2],type="l",main="Death",xlab="Patient No.",
     ylab="CUSUM Statistic")
    lines(c(0,1e99), c(70,70),lty=2,col=2)
    lines(c(0,1e99), c(38,38),lty=2,col=3)

######################################################################
# Run the CUSUM with thresholds as in Steiner et al. (1999).
# After each alarm the CUSUM statistic is set to zero and
# monitoring continues from this point. Triangles indicate alarm
# in the respective CUSUM (nearmiss or death). If in both
# simultaneously then an alarm is caued by the secondary limits.
######################################################################
pb.surv2 <- pairedbinCUSUM(deleval,control=list(theta0=theta0,
             theta1=theta1,h1=70*w1,h2=32*w2,h11=38*w1,h22=17*w2))

plot(pb.surv2, xaxis.labelFormat=NULL)
}

\author{S. Steiner and M. \enc{Höhle}{Hoehle}}
\keyword{regression}
