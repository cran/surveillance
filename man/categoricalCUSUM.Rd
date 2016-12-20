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

  Once an alarm is found the CUSUM scheme is resetted (to zero) and
  monitoring continues from there.
}
\seealso{\code{\link{categoricalCUSUM}}}
\value{An \code{sts} object with \code{observed}, \code{alarm},
  etc. slots trimmed to the \code{control$range} indices.
}
\references{
  Höhle, M. (2010):
  Online Change-Point Detection in Categorical Time Series.
  Book chapter in T. Kneib and G. Tutz (Eds.), Statistical
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

  abattoir2 <- new("sts",epoch=1:nrow(abattoir), start=c(2006,1),freq=52,
    observed=cbind(abattoir@observed,abattoir@populationFrac -abattoir@observed),
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

#ToDo: Compute run length using LRCUSUM.runlength
}

\author{M. \enc{Höhle}{Hoehle}}
\keyword{regression}

% positives <- matrix(c(25,54,50,70,54,83,62,37,53,29,48,95,63,53,47,67,31,56,20,37,33,34,28,58,45,31,34,17,25,23,19,32,29,34,46,58,45,50,50,41,57,66,41,68,58,35,40,45,20,49,33,24,21,31,25,28,29,20,25,13,41,31,25,25,7,9,48,15,4,45,7,52,21,19,60,25,24,26,16,15,8,17,43,19,10,31,26,19,20,20,27,17,27,23,35,34,53,32,58,81,37,3,0,0,0,1,41,21,34,33,48,17,61,39,70,13,56,65,50,39,46,29,35,18,35,43,34,49,49,18,36,30,34,55,61,44,35,23,35,37,33,37,15,17,66,8,2,10,39,17,42,27,28,55,26,0,0,0,0,0,1,0,0,0,16,36,36,25,34,26,2),ncol=1)

% n <- matrix(c(4200,4900,4700,5100,4900,5000,4700,4500,4500,4400,4300,4500,4500,4800,3000,3600,4200,4300,2900,4000,3000,4500,3300,4700,4300,4100,3100,2500,3700,3300,3700,4200,4200,4100,4000,5000,3900,3700,4600,4400,4500,5000,4100,3800,2300,3000,3700,3200,2000,2500,4000,1500,1700,3000,3000,2500,2300,3100,2600,1300,4100,1800,2600,2900,1100,1300,3700,2700,1700,2700,1100,3300,1800,3500,4100,3300,4100,2500,3000,2400,800,3300,3100,3000,3100,2500,3300,3200,2700,2700,1500,1400,2200,1600,2700,1700,4200,3900,8100,5900,4800,300,0,0,0,200,5000,5400,6100,5500,4300,4100,8000,8200,7700,2500,5200,7200,7400,5000,6500,3900,5300,4100,6600,6600,5400,6300,6100,5000,5000,4700,4900,5200,4800,4900,4900,4900,5000,5000,4900,4900,900,800,1300,200,700,1200,2300,1600,2900,2600,3600,4000,2000,0,0,0,0,0,100,0,0,100,2200,2900,3100,3200,3900,4300,1000),ncol=1)

%#Create S4 object
%sts <- new("sts",epoch=1:nrow(positives),start=c(2006,1),freq=52,observed=positives,state=0*positives,populationFrac=n,multinomialTS=TRUE)
