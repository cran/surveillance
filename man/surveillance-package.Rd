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
  categorical data. Focus is on outbreak detection in count data time
  series originating from public health surveillance of infectious
  diseases, but applications could just as well originate from
  environmetrics, reliability engineering, econometrics or social
  sciences.
 
} 

\details{
  \tabular{ll}{
    Package: \tab surveillance\cr
    Type: \tab Package\cr
    Version: \tab 1.1-0\cr
    Date: \tab 2009-10-14\cr
    License: \tab GPL version 2 (http://www.gnu.org/licenses/gpl.html) \cr
  }

  \code{surveillance} is an \R package implementing statistical methods for the
  retrospective modeling and prospective change-point detection in time
  series of counts, proportions and categorical data. The main
  application is in the detection of aberrations in routine collected
  public health data seen as univariate and multivariate time series of
  counts, but applications could just as well originate from
  environmetrics, econometrics or social sciences. As many methods rely
  on statistical process control methodology, the package is thus also
  relevant to quality control and reliability engineering.

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
 model in Held et al. (2005) and Paul et al. (2008) handling
 multivariate time series of counts. Finally, the fully Bayesian
 approach for univariate time series of counts from Held et al. (2006)
 is also implemented.
}
\author{
  Author: M. \enc{Höhle}{Hoehle} with contributions from T. Correa, M. Hofmann,
  C. Lang, M. Paul, A. Riebler, S. Steiner and V. Wimmer
  
  Maintainer: Michael \enc{Höhle}{Hoehle} \email{<hoehle@stat.uni-muenchen.de>}
}
\references{
  surveillance: An R package for the surveillance of infectious diseases (2007), M. \enc{Höhle}{Hoehle}, Computational Statistics, 22(4), pp. 571---582.
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
