\encoding{latin1}
\name{surveillance-package}
\alias{surveillance-package}
\alias{surveillance}
\docType{package}

\title{
  Temporal and Spatio-Temporal Modeling and Monitoring of Epidemic Phenomena
%  Outbreak Detection and Retrospective Modeling for Surveillance Data
}

\description{
  A package implementing statistical methods for the modeling and
  change-point detection in time series of counts, proportions and
  categorical data, as well as for the modeling of continuous-time
  epidemic phenomena, e.g. discrete-space setups such as the spatially
  enriched Susceptible-Exposed-Infectious-Recovered (SEIR) models for
  surveillance data, or continuous-space point process data such as the
  occurrence of disease or earthquakes. Main focus is on outbreak
  detection in count data time series originating from public health
  surveillance of infectious diseases, but applications could just as well
  originate from environmetrics, reliability engineering, econometrics or
  social sciences.
  Currently the package contains implementations of typical outbreak
  detection procedures such as Stroup et. al (1989), Farrington et al
  (1996), Rossi et al (1999), Rogerson and Yamada (2001), a Bayesian
  approach, negative binomial CUSUM methods and a detector based on
  generalized likelihood ratios.
  Furthermore, inference methods for the retrospective infectious disease
  model in Held et al (2005), Held et al (2006), Paul et al (2008) and
  Paul and Held (2011) are provided. A novel CUSUM approach combining
  logistic and multinomial logistic modelling is also included.
  Continuous self-exciting spatio-temporal point processes are modeled
  through additive-multiplicative conditional intensities as described in
  H\enc{ö}{oe}hle (2009) ("twinSIR", discrete space) and Meyer et al (2012)
  ("twinstim", continuous space).
  The package contains several real-world datasets, the ability to
  simulate outbreak data, visualize the results of the monitoring in
  temporal, spatial or spatio-temporal fashion.
} 

\details{
  \tabular{ll}{
    Package: \tab surveillance\cr
    Version: \tab 1.4-3\cr
    Date: \tab 2012-10-01\cr
    License: \tab GPL version 2
    (\url{http://www.gnu.org/licenses/old-licenses/gpl-2.0.html}) \cr
    URL: \tab \url{http://surveillance.r-forge.r-project.org/}\cr
  }
  
  \pkg{surveillance} is an \R package implementing statistical methods
  for the retrospective modeling and prospective change-point detection
  in time series of counts, proportions and categorical data. The main
  application is in the detection of aberrations in routine collected
  public health data seen as univariate and multivariate time series of
  counts or point-processes. However, applications could just as well
  originate from environmetrics, econometrics or social sciences. As
  many methods rely on statistical process control methodology, the
  package is thus also relevant to quality control and reliability
  engineering.

  The fundamental data structure of the package is an S4 class
  \code{sts} wrapping observations, monitoring results and date handling
  for multivariate time series. Currently the package contains
  implementations typical outbreak detection procedures such as Stroup et
  al. (1989), Farrington et al., (1996), Rossi et al. (1999), Rogerson
  and Yamada (2001), a Bayesian approach (H\enc{ö}{oe}hle, 2007),
  negative binomial CUSUM methods (H\enc{ö}{oe}hle and Mazick, 2009), and a
  detector based on generalized likelihood ratios (H\enc{ö}{oe}hle
  and Paul, 2008). However, also CUSUMs for the prospective change-point
  detection in binomial, beta-binomial and multinomial time series is
  covered based on generalized linear modelling.  This includes
  e.g. paired binary CUSUM described by Steiner et al. (1999) or paired
  comparison Bradley-Terry modelling described in H\enc{ö}{oe}hle
  (2010). The package contains several real-world datasets, the ability
  to simulate outbreak data, visualize the results of the monitoring in
  temporal, spatial or spatio-temporal fashion.

  Furthermore, inference methods for the retrospective infectious disease
  model in Held et al. (2005), Paul et al. (2008) ("algo.hhh") and Paul
  and Held (2011) ("hhh4") handling multivariate time series of
  counts. Furthermore, the fully Bayesian approach for univariate time
  series of counts from Held et al. (2006) ("twins") is also implemented
  Self-exciting spatio-temporal point processes are modeled through
  additive-multiplicative conditional intensities as described in
  H\enc{ö}{oe}hle (2009) ("twinSIR") and Meyer et al (2012) ("twinstim").

  Altogether, the package contains several real-world datasets, the
  ability to simulate outbreak data, visualize the results of the
  monitoring in temporal, spatial or spatio-temporal fashion.
}

\section{Acknowledgements}{
  T. Correa, L. Held, M. Hofmann, C. Lang, A. Riebler,
  D. Saban\enc{é}{e}s Bov\enc{é}{e}, S. Steiner, M. Virtanen,
  and V. Wimmer made substantial contributions of code.

  The authors would like to thank the following people for ideas,
  discussions, testing and feedback:
  Doris Altmann, Johannes Dreesman, Johannes Elias, Mayeul Kauffmann,
  Marc Geilhufe, Marcos Prates, Barry Rowlingson, Christopher W. Ryan,
  Klaus Stark, Yann Le Strat, Andr\enc{é}{e} Michael Toschke,
  Achim Zeileis.
}

\author{
  Michael H\enc{ö}{oe}hle, Sebastian Meyer, and Michaela Paul
  
  Maintainer: Michael H\enc{ö}{oe}hle
  \email{Michael.Hoehle@stat.uni-muenchen.de}
}

\references{
  See \code{citation(package="surveillance")}.
}

\keyword{ package }

\examples{
#Code from an early survey article about the package: Hoehle (2007)
#available from http://surveillance.r-forge.r-project.org/
\dontrun{demo(cost)}
#Code from a more recent book chapter about using the package for the
#monitoring of Danish mortality data (Hoehle, 2009).
\dontrun{demo(biosurvbook)}
}
