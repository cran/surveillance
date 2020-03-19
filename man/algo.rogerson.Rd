\name{algo.rogerson}
\alias{algo.rogerson}
\title{Modified CUSUM method as proposed by Rogerson and Yamada (2004)}
\description{
  Modified Poisson CUSUM method that allows for a time-varying in-control parameter
  \eqn{\theta_{0,t}} as proposed by Rogerson and Yamada (2004). The
  same approach can be applied to binomial data if \code{distribution="binomial"}
  is specified.
}
\usage{
algo.rogerson(disProgObj, control = list(range = range,
   theta0t = NULL, ARL0 = NULL, s = NULL, hValues = NULL,
   distribution = c("poisson","binomial"), nt = NULL, FIR=FALSE, 
   limit = NULL, digits = 1))
}
\arguments{
  \item{disProgObj}{object of class \code{disProg} that includes a matrix with
     the observed number of counts}
  \item{control}{ list with elements
    \describe{
      \item{range}{vector of indices in the observed matrix of \code{disProgObj} to monitor}
      \item{theta0t}{matrix with in-control parameter, must be specified}
      \item{ARL0 }{ desired average run length \eqn{\gamma} }
      \item{s}{change to detect, see \code{\link{findH}} for further details}
      \item{hValues}{matrix with decision intervals \code{h} for a sequence of
           values \eqn{\theta_{0,t}}
           (in the range of \code{theta0t}) }
      \item{distribution}{\code{"poisson"} or \code{"binomial"} }
      \item{nt}{optional matrix with varying sample sizes for the binomial CUSUM}
      \item{FIR}{a FIR CUSUM with head start \eqn{\frac{\code{h}}{2}}
            is applied to the data if \code{TRUE}, otherwise no head start is used;
            see details }
      \item{limit}{numeric that determines the procedure after an alarm is given, see details}
      \item{digits}{the reference value and decision interval are rounded to
         \code{digits} decimal places. Defaults to 1 and should correspond
         to the number of digits used to compute \code{hValues} }
    }
  }
}
\details{
The CUSUM for a sequence of Poisson or binomial
variates \eqn{x_t} is computed as
\deqn{S_t = \max \{0, S_{t-1} + c_t (x_t- k_t)\} , \, t=1,2,\ldots ,}
where \eqn{S_0=0} and \eqn{c_t=\frac{h}{h_t} }; \eqn{k_t} and \eqn{h_t} are time-varying
reference values and decision intervals.
An alarm is given at time \eqn{t} if \eqn{S_t \geq h}.

If \code{FIR=TRUE}, the CUSUM starts
with a head start value \eqn{S_0=\frac{\code{h}}{2}} at time \eqn{t=0}.
After an alarm is given, the FIR CUSUM starts again at this head start value.
  
The procedure after the CUSUM gives an alarm can be determined by \code{limit}.
Suppose that the CUSUM signals at time \eqn{t}, i.e. \eqn{S_t \geq h}.
For numeric values of \code{limit}, the CUSUM is bounded
above after an alarm is given, % at time \eqn{t-1},
i.e. \eqn{S_{t}} is set to \eqn{ \min\{\code{limit} \cdot h,S_{t}\} }.
%\deqn{S_{t} = \max \{0, S_{t-1} + c_t(x_t - k_t)\}. }
Using \code{limit}=0 corresponds to
resetting \eqn{S_t} to zero after an alarm as proposed in the original
formulation of the CUSUM. If \code{FIR=TRUE},
\eqn{S_{t}} is reset to \eqn{ \frac{\code{h}}{2} }
(i.e. \code{limit}=\eqn{\frac{\code{h}}{2} } ).
If \code{limit=NULL}, no resetting occurs after an alarm is given.
}
\note{\code{algo.rogerson} is a univariate CUSUM method. If the data are
available in several regions (i.e. \code{observed} is a matrix),
multiple univariate CUSUMs are applied to each region.
}
\value{Returns an object of class \code{survRes} with elements
 \item{alarm}{indicates whether the CUSUM signaled at time \eqn{t} or not
            (1 = alarm, 0 = no alarm) }
 \item{upperbound}{CUSUM values \eqn{S_{t}} }
 \item{disProgObj}{\code{disProg} object }
 \item{control}{list with the alarm threshold \eqn{h} and the specified
    control object}
}
\examples{
# simulate data (seasonal Poisson)
set.seed(123)
t <- 1:300
lambda <- exp(-0.5 + 0.4 * sin(2*pi*t/52) + 0.6 * cos(2*pi*t/52))
data <- create.disProg(week = t, observed = rpois(length(lambda), lambda))

# determine a matrix with h values
hVals <- hValues(theta0 = 10:150/100, ARL0=500, s = 1, distr = "poisson")

# apply modified Poisson CUSUM
res <- algo.rogerson(data, control=c(hVals, list(theta0t=lambda,range=1:300)))
plot(res)
}
\references{
Rogerson, P. A. and Yamada, I. Approaches to Syndromic Surveillance When Data
Consist of Small Regional Counts. Morbidity and Mortality Weekly Report, 2004,
53/Supplement, 79-85
}
\seealso{\code{\link{hValues}}}
\keyword{classif}
