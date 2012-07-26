\encoding{latin1}
\name{surveillance-package}
\alias{surveillance-package}
\alias{surveillance}
\docType{package}
\title{
  Outbreak detection algorithms for surveillance data
}
\description{

  A package implementing statistical methods for the modeling and
  change-point detection in time series of counts, proportions and
  categorical data. Focus is on the modelling and change-point detection
  in count data time series or point-processes originating from public
  health surveillance of infectious diseases, but applications could
  just as well originate from environmetrics, reliability engineering,
  econometrics or social sciences.
 
} 

\details{
  \tabular{ll}{
    Package: \tab surveillance\cr
    Type: \tab Package\cr
    Version: \tab 1.1-0\cr
    Date: \tab 2009-10-14\cr
    License: \tab GPL version 2 (http://www.gnu.org/licenses/gpl.html) \cr
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
 and Yamada (2001), a Bayesian approach (\enc{Höhle}{Hoehle}, 2007),
 negative binomial CUSUM methods (\enc{Höhle}{Hoehle} and Mazick, 2009), and a
 detector based on generalized likelihood ratios (\enc{Höhle}{Hoehle}
 and Paul, 2008). However, also CUSUMs for the prospective change-point
 detection in binomial, beta-binomial and multinomial time series is
 covered based on generalized linear modelling.  This includes
 e.g. paired binary CUSUM described by Steiner et al. (1999) or paired
 comparison Bradley-Terry modelling described in \enc{Höhle}{Hoehle}
 (2010). The package contains several real-world datasets, the ability
 to simulate outbreak data, visualize the results of the monitoring in
 temporal, spatial or spatio-temporal fashion.

 Furthermore, inference methods for the retrospective infectious disease
 model in Held et al. (2005), Paul et al. (2008) ("algo.hhh") and Paul
 and Held (2011) ("hhh4") handling multivariate time series of
 counts. Furthermore, the fully Bayesian approach for univariate time
 series of counts from Held et al. (2006) ("twins") is also implemented
 Self-exciting spatio-temporal point processes are modeled through
 additive-multiplicative conditional intensities as described in Höhle
 (2009) ("twinSIR") and Meyer et al (2012) ("twinstim").

 Altogether, the package contains several real-world datasets, the
 ability to simulate outbreak data, visualize the results of the
 monitoring in temporal, spatial or spatio-temporal fashion.

}
\author{
  Author: M. H\enc{ö}{oe}hle, S. Meyer and M. Paul
  with contributions from T. Correa, M. Hofmann, C. Lang, A. Riebler,
  D. Sabanés Bové, S. Steiner, M. Virtanen and V. Wimmer
  
  Maintainer: Michael H\enc{ö}{oe}hle \email{<hoehle@stat.uni-muenchen.de>}
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
