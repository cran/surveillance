\name{categoricalCUSUM}
\alias{categoricalCUSUM}
\alias{catcusum.LLRcompute}
\encoding{latin1}

\title{CUSUM detector for time-varying categorical time series}
\description{
  Function to process \code{sts} object by binomial, beta-binomial
  or multinomial CUSUM as described by \enc{Höhle}{Hoehle} (2010).
  Logistic, multinomial logistic, proportional
  odds or Bradley-Terry regression models are used to specify in-control
  and out-of-control parameters.
  The implementation is illustrated in Salmon et al. (2016).
}
\usage{
categoricalCUSUM(stsObj,control = list(range=NULL,h=5,pi0=NULL,
                 pi1=NULL, dfun=NULL, ret=c("cases","value")),...)
}
\arguments{
  \item{stsObj}{Object of class \code{sts} containing the number of
    counts in each of the \eqn{k} categories of the response
    variable. Time varying number of counts \eqn{n_t} is found in slot
    \code{populationFrac}. }
  \item{control}{Control object containing several items
    \itemize{
      \item{\code{range}}{Vector of length \eqn{t_{max}} with indices of the
	\code{observed} slot to monitor.}
      \item{\code{h}}{Threshold to use for the monitoring. Once the
	CUSUM statistics is larger or equal to \code{h} we have an alarm.}
      \item{\code{pi0}}{\eqn{(k-1) \times t_{max}} in-control probability
	vector for all categories except the reference category.}
      \item{\code{mu1}}{\eqn{(k-1) \times t_{max}} out-of-control probability
	vector for all categories except the reference category.}
      \item{\code{dfun}}{The probability mass function (PMF) or density used
	to compute the likelihood ratios of the CUSUM. In a negative
	binomial CUSUM this is \code{dnbinom}, in a binomial CUSUM
	\code{dbinom} and in a multinomial CUSUM \code{dmultinom}. The
	function must be able to handle the arguments \code{y},
	\code{size}, \code{mu} and \code{log}. As a consequence, one in
	the case of, e.g, the beta-binomial distribution has to write a small
	wrapper function.}
      \item{\code{ret}}{Return the necessary proportion to sound an alarm in the
	slot \code{upperbound} or just the value of the CUSUM
	statistic. Thus, \code{ret} is one of the values in
	\code{c("cases","value")}. Note: For the binomial PMF it is
	possible to compute this value explicitly, which is much faster
	than the numeric search otherwise conducted. In case \code{dfun}
	just corresponds to \code{dbinom} just set the attribute
	\code{isBinomialPMF} for the \code{dfun} object.}
  }}
  \item{\dots}{Additional arguments to send to \code{dfun}.}
}
\details{
  The function allows the monitoring of categorical time series as
  described by regression models for binomial, beta-binomial or
  multinomial data. The later includes e.g. multinomial logistic
  regression models, proportional odds models or Bradley-Terry models
  for paired comparisons. See the \enc{Höhle}{Hoehle} (2010) reference
  for further details about the methodology.

  Once an alarm is found the CUSUM scheme is reset (to zero) and
  monitoring continues from there.
}
\seealso{\code{\link{LRCUSUM.runlength}}}
\value{An \code{sts} object with \code{observed}, \code{alarm},
  etc. slots trimmed to the \code{control$range} indices.
}
\references{
  \enc{Höhle}{Hoehle}, M. (2010):
  Online Change-Point Detection in Categorical Time Series.
  In: T. Kneib and G. Tutz (Eds.), Statistical
  Modelling and Regression Structures, Physica-Verlag.

  Salmon, M., Schumacher, D. and \enc{Höhle}{Hoehle}, M. (2016):
  Monitoring count time series in \R: Aberration detection in public
  health surveillance. \emph{Journal of Statistical Software},
  \bold{70} (10), 1-35. \doi{10.18637/jss.v070.i10}
}
\examples{
if (require("gamlss")) {
  ###########################################################################
  #Beta-binomial CUSUM for a small example containing the time-varying
  #number of positive test out of a time-varying number of total
  #test.
  #######################################

  #Load meat inspection data
  data("abattoir")

  #Use GAMLSS to fit beta-bin regression model
  phase1 <- 1:(2*52)
  phase2  <- (max(phase1)+1) : nrow(abattoir)

  #Fit beta-binomial model using GAMLSS
  abattoir.df <- as.data.frame(abattoir)

  #Replace the observed and epoch column names to something more convenient
  dict <- c("observed"="y", "epoch"="t", "population"="n")
  replace <- dict[colnames(abattoir.df)]
  colnames(abattoir.df)[!is.na(replace)] <- replace[!is.na(replace)]

  m.bbin <- gamlss( cbind(y,n-y) ~ 1 + t +
                      + sin(2*pi/52*t) + cos(2*pi/52*t) +
                      + sin(4*pi/52*t) + cos(4*pi/52*t), sigma.formula=~1,
                   family=BB(sigma.link="log"),
                   data=abattoir.df[phase1,c("n","y","t")])

  #CUSUM parameters
  R <- 2 #detect a doubling of the odds for a test being positive
  h <- 4 #threshold of the cusum

  #Compute in-control and out of control mean
  pi0 <- predict(m.bbin,newdata=abattoir.df[phase2,c("n","y","t")],type="response")
  pi1 <- plogis(qlogis(pi0)+log(R))
  #Create matrix with in control and out of control proportions.
  #Categories are D=1 and D=0, where the latter is the reference category
  pi0m <- rbind(pi0, 1-pi0)
  pi1m <- rbind(pi1, 1-pi1)


  ######################################################################
  # Use the multinomial surveillance function. To this end it is necessary
  # to create a new abattoir object containing counts and proportion for
  # each of the k=2 categories. For binomial data this appears a bit
  # redundant, but generalizes easier to k>2 categories.
  ######################################################################

  abattoir2 <- sts(epoch=1:nrow(abattoir), start=c(2006,1), freq=52,
    observed=cbind(abattoir@observed, abattoir@populationFrac-abattoir@observed),
    populationFrac=cbind(abattoir@populationFrac,abattoir@populationFrac),
    state=matrix(0,nrow=nrow(abattoir),ncol=2),
    multinomialTS=TRUE)

  ######################################################################
  #Function to use as dfun in the categoricalCUSUM
  #(just a wrapper to the dBB function). Note that from v 3.0-1 the
  #first argument of dBB changed its name from "y" to "x"!
  ######################################################################
  mydBB.cusum <- function(y, mu, sigma, size, log = FALSE) {
    return(dBB(y[1,], mu = mu[1,], sigma = sigma, bd = size, log = log))
  }


  #Create control object for multinom cusum and use the categoricalCUSUM
  #method
  control <- list(range=phase2,h=h,pi0=pi0m, pi1=pi1m, ret="cases",
		   dfun=mydBB.cusum)
  surv <- categoricalCUSUM(abattoir2, control=control,
			   sigma=exp(m.bbin$sigma.coef))

  #Show results
  plot(surv[,1],dx.upperbound=0)
  lines(pi0,col="green")
  lines(pi1,col="red")

  #Index of the alarm
  which.max(alarms(surv[,1]))
}
}

\author{M. \enc{Höhle}{Hoehle}}
\keyword{regression}
