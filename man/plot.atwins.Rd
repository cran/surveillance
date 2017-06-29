\name{plot.atwins}
\alias{plot.atwins}
\encoding{latin1}
\title{Plot results of a twins model fit} 
\description{
  Plot results of fitting a twins model using MCMC output. Plots similar
  to those in the Held et al. (2006) paper are generated
}
\usage{
  \method{plot}{atwins}(x, which=c(1,4,6,7), ask=TRUE, \dots) 
}
\arguments{
  \item{x}{An object of class \code{atwins}.}
  \item{which}{a vector containing the different plot types to show
    \describe{
      \item{1}{A plot of the observed time series Z is shown together
	with posterior means for the number of endemic cases (X) and number of
	epidemic cases (Y).}
      \item{2}{This plot shows trace plots of the gamma parameters
	over all MCMC samples.}
      \item{3}{This shows a trace plot of psi, which controls the
	overdispersion in the model.}
      \item{4}{Autocorrelation functions for K and psi are shown in
	order to judge whether the MCMC sampler has converged.}
      \item{5}{Shows a plot of the posterior mean of the seasonal
	model nu[t] together with 95\% credibility intervals based on the
	quantiles of the posterior.}
      \item{6}{Histograms illustrating the posterior density for K and
	psi. The first one corresponds to Fig. 4(f) in the paper.}
      \item{7}{Histograms illustrating the predictive posterior
	density for the next observed number of cases Z[n+1]. Compare with Fig.5
	in the paper.}
    }
  }
  \item{ask}{Boolean indicating whether to ask for a newline before
    showing the next plot.}
  \item{\dots}{Additional control for the plots, which are currently ignored.}
}
\details{
  For details see the plots in the paper. Basically MCMC output is
  visualized. This function is together with \code{algo.twins} still
  experimental.  
}
\value{This function does not return anything.}
\references{Held, L., Hofmann, M., \enc{Höhle}{Hoehle}, M. and Schmid V. (2006)
  A two-component
  model for counts of infectious diseases, Biostatistics, \bold{7}, pp.
  422--437.
}
\author{M. Hofmann and M. \enc{Höhle}{Hoehle}}
\seealso{\link{algo.twins}}
\examples{
\dontrun{
#Apparently, the algo.atwins can crash on some LINUX systems
#thus for now the example section is commented

#Load the data used in the Held et al. (2006) paper
data("hepatitisA")

#Fix seed - this is used for the MCMC samplers in twins
set.seed(123)

#Call algorithm and save result
otwins <- algo.twins(hepatitisA)

#This shows the entire output
plot(otwins,which=c(1,2),ask=FALSE)
}
}
\keyword{ts}
\keyword{regression}
