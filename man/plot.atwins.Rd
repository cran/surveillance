\name{plot.atwins}
\alias{plot.atwins}
\encoding{latin1}
\title{Plots for Fitted \code{algo.twins} Models}
\description{
  Plot results of fitting a twins model using MCMC output. Plots similar
  to those in the Held et al. (2006) paper are generated.
}
\usage{
  \method{plot}{atwins}(x, which=c(1,4,6,7), ask=TRUE, \dots) 
}
\arguments{
  \item{x}{An object of class \code{"atwins"} as returned by
    \code{\link{algo.twins}}.}
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
    showing the next plot (only if multiple are shown).}
  \item{\dots}{Additional arguments for \code{\link{stsplot_time}},
    used for plot type 1.}
}
\details{
  For details see the plots in the paper. Basically MCMC output is
  visualized. This function is experimental, as is \code{\link{algo.twins}}.
}
%\value{NULL (invisibly)}
\references{Held, L., Hofmann, M., \enc{Höhle}{Hoehle}, M. and Schmid V. (2006)
  A two-component
  model for counts of infectious diseases, Biostatistics, \bold{7}, pp.
  422--437.
}
\author{M. Hofmann and M. \enc{Höhle}{Hoehle}}
\seealso{\code{\link{algo.twins}} (with an example)}
\keyword{hplot}
